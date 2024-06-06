% In the gp2 model analysis, the models are fit to the observed data points
% and the responses. In contrast, gp1 fits the models to the responses only
%
% The idea is that the gp simultaneously estimates parameters for the
% presented graph but the fit takes into account the predicted probability
% of the response

clear all
clc
close all

subjects = dlmread('socfun_subject_numbers.txt');

% matlabpool local 4
totalTime = 0; 
for sidx2 = 1 % numel(subjects)
    keep sidx2 subjects totalTime
    load(fullfile(pwd, 'Data', sprintf('final_%d_data.mat', subjects(sidx2))))
    data = finalDrawData;
    data(ismember(data(:,2), [9 21]), :) = []; % Remove sine functions
        
    tic
    
    startfrom = 1;
    startend = 1; %24 ;
    trialidx = [1:8 10:20 22:26]';
    
    logmodelprior = repmat(log(1/5), 1, 5);
    
    for i = startfrom:startend
        disp([sidx2 i])
        
        %% Set up data for fitting
        % Observed data
        X     = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,4)), 4); % Observed X
        Y     = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,5)), 5); % Observed Y
        
        % Response
        XSTARall = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,6)), 6); % Response X
        YSTARall = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,7)), 7); % Response Y
        
        YSTARall(XSTARall < min(X) | XSTARall > max(X)) = []; 
        XSTARall(XSTARall < min(X) | XSTARall > max(X)) = [];
        
        [XSTAR, idx, jdx] = unique(XSTARall);
        YSTAR = YSTARall(idx,1);
        
        % Reorder data in case something is misaligned
        [xy, xyord] = sortrows([X, Y], 1);
        xystar      = sortrows([XSTAR YSTAR], 1);
        
        % Observed data put in order of X
        X = xy(:,1);
        Y = xy(:,2);
        
        % Response data put in order of X
        XSTAR = xystar(:,1);
        YSTAR = xystar(:,2);
       
        %% GP analysis for model 1
        if ~isempty(X) && ~isempty(XSTAR)
            % Sample 40 responses evenly spaced along the length of the response
%             respLocations = unique(floor(linspace(1, numel(XSTAR), 40)));
            rlsets = linspace(min(X), max(X), 40);
            for rlIdx = 1:numel(rlsets)
                respLocations(rlIdx) = find(abs(rlsets(rlIdx) - XSTAR) == min(abs(rlsets(rlIdx) - XSTAR)), 1, 'first');
            end

            xSTAR = XSTAR(respLocations);
            ySTAR = YSTAR(respLocations);
            
            maxi = 5;
            for modidx = 1:maxi
                if modidx  <= (maxi - 1)
                    k = modidx ; % Degree of polynomial to fit
                    covfunc = {'covSum', {sprintf('covPoly%d', k),  'covNoise'}}; % Set up covariance function - polynomial basis function of degree k with additive noise

                    logtheta = [-ones(1, k+1), -1]';
                    parms{i,modidx} = minimize(logtheta, 'gpr2', -100, covfunc, X, Y, xSTAR, ySTAR, 2);
                    
                    [lnml dnml] = gpr2(parms{i,modidx}, covfunc, X, Y, xSTAR, ySTAR, 2);      % Compute the fit of the model - p(y|M) (marginal including hyperparameters)
                    [fmu FS2] = gpr2(parms{i,modidx}, covfunc, X, Y, xSTAR, ySTAR, 1); % Compute the mean function and standard deviation for prediction space
                    fS2 = FS2 - exp(2 * parms{i,modidx}(end)); % Subtract off additive noise
                    
                    M = numel(logtheta); % eval(feval(covfunc{:}));
                    eparms = [exp(-2 * parms{i, modidx}(1:(M - 1))); exp(2 * parms{i, modidx}(M))]';
                    lnPrior = logmvnpdf(eparms, zeros(1,M), diag(repmat(100, 1, M)));
                    [lnml, dnml, H] = gpr2(parms{i,modidx}, covfunc, X, Y, xSTAR, ySTAR, 0);
                else
                    covfunc = {'covSum', {'covSEiso', 'covNoise'}}; % Set up covariance function - polynomial basis function of degree k with additive noise
                    
                    logtheta = [-1 -1 -1]';
                    parms{i,modidx} = minimize(logtheta, 'gpr2', -100, covfunc, X, Y, xSTAR, ySTAR, 2);  % Maximize the marginal probability of the data with respect to the hyperparamters
                    
                    [lnml dnml] = gpr2(parms{i, modidx}, covfunc, X, Y, xSTAR, ySTAR, 2);      % Compute the fit of the model - p(y|M) (marginal including hyperparameters)
                    [fmu FS2] = gpr2(parms{i, modidx}, covfunc, X, Y, xSTAR, ySTAR, 1); % Computer the mean function and standard deviation for prediction space
                    fS2 = FS2 - exp(2 * parms{i, modidx}(end)); % Subtract off additive noise
                    
                    M = numel(logtheta); % eval(feval(covfunc{:}));
                    eparms = [exp(parms{i, modidx}(1)); exp(2 * parms{i, modidx}(2)); exp(2 * parms{i, modidx}(3))]';
                    lnPrior = logmvnpdf(eparms, zeros(1,M), diag(repmat(100, 1, M)));
                    [lnml, dnml, H] = gpr2(parms{i,modidx}, covfunc, X, Y, xSTAR, ySTAR, 0);
                end
                lpD(modidx) = -lnml - .5 * log(abs(det(H))) + M/2 * log(2 * pi) + lnPrior;
                
                K = feval(covfunc{:}, parms{i, modidx}, X);
                % [E, v] = eigDecomp(K, 1/exp(2 * floghyper(end)), -5, 5);
                
                means{modidx} = fmu;
                vars{modidx} = fS2;
                nParms(modidx) = numel(parms{i, modidx});
                fits(modidx) = -lnml;
            end
            fullfits.means(i,:)  = means;
            fullfits.vars(i,:)   = vars;
            fullfits.lpD(i,:)    = lpD;
            fullfits.nParms(i,:) = nParms;
            fullfits.fits(i,:)   = fits;
            fullfits.pM(i,:)     = (exp(fullfits.lpD(i, :)) .* exp(logmodelprior))/(exp(fullfits.lpD(i,:)) * exp(logmodelprior)');
        end
        fullfits.X{i} = X;
        fullfits.Y{i} = Y;
        fullfits.xSTAR{i} = xSTAR;
        fullfits.ySTAR{i} = ySTAR;
    end
    fullfits.parms       = parms;
    t = toc; 
    displayTime(t);
    totalTime = totalTime + t;
    save(fullfile(pwd, sprintf('gpAnalysis_Full_s%d.mat', subjects(sidx2))))
%     clear fullfits
end
%  matlabpool close
displayTime(totalTime);
 % figure; fill([fullfits.xSTAR{i}; flipud(fullfits.xSTAR{i})], [fullfits.means{i,j} - fullfits.vars{i,j}; flipud(fullfits.means{i,j}) + fullfits.vars{i,j}], [.75 .75 .75]); hold on;plot(fullfits.xSTAR{i}, fullfits.ySTAR{i}, ' ok', fullfits.xSTAR{i}, fullfits.means{i,j}, '-k', fullfits.X{i}, fullfits.Y{i}, ' *k')
  