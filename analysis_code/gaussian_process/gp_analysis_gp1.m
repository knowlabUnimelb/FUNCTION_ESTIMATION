% In the gp2 model analysis, the models are fit to the observed data points
% and the responses. In contrast, gp1 fits the models to the responses only
%
% This could potentially fail because the variance in the responses is low.
% For instance, a hand-drawn straight line may not be perfectly straight
% but may have some curvature.  This curvature is problematic and would be
% ignored if the drawn line has any variability at all.  With no
% variability, a more complex model than necessarity is likely to fit
% better (even with penalties applied for bad fits).
%
% To get around this, I'm fitting a single model (squared exponential
% kernel) and using the parameters of that kernel to differentiate
% polyonmials from data tracking
%
% uses gpml-toolbox-v1
 
clear all
clc
close all
reset(symengine)

subjects = dlmread('socfun_subject_numbers.txt');
subjects(ismember(subjects,  [75 101])) = [];
% subjects = (171:181)';

% matlabpool local 4

for sidx2 = 1%:numel(subjects) % Cycle through each of the subjects
    keep sidx2 subjects % Everything is saved to a structure each loop so delete these variables
    
    % Load the data
    load(fullfile(pwd, 'Data', sprintf('final_%d_data.mat', subjects(sidx2))))
    data = finalDrawData;
    data(ismember(data(:,2), [9 21]), :) = []; % Remove sine functions
    
    tic
    
    % Select trials to examine. 
    startfrom = 1;
    startend =  24;
    trialidx = [1:8 10:20 22:26]'; % Trials 9 and 21 are removed as sine functions
    
    % Models
    models = {'covPoly1', 'covPoly2', 'covPoly3', 'covPoly4', 'covSEiso'};
    
    % Other parameters
    nRespLocations = 40; % Number of points to use for each response
    logmodelprior = repmat(log(1/5), 1, 5); % Prior over models (i.e., uniform)
    
    % Cycle through functions for this subject
    for i = startfrom:startend
        disp([sidx2 i])
        
        %% Set up data for fitting
        
        % Observed data
        X     = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,4)), 4); % Observed X
        Y     = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,5)), 5); % Observed Y
        
        % Response
        XSTARall = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,6)), 6); % Response X
        YSTARall = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,7)), 7); % Response Y
        
        % Select out unique X points only
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
        
        
        %% GP analysis
        % If the current data set or response is not empty (which is should not be if a response exists for this item), then fit the models
        if ~isempty(X) && ~isempty(XSTAR) && numel(XSTAR) > 1
            
            % Sample 40 responses evenly spaced along the length of the response
            respLocations = unique(floor(linspace(1, numel(XSTAR), nRespLocations)));
            xSTAR = XSTAR(respLocations);
%             ySTAR = YSTAR(respLocations) + randn(nRespLocations, 1)/10;
            ySTAR = YSTAR(respLocations);
            
            % Get predictions for these x values
            xPRED = linspace(-1,1,50)';
            
            % Cycle through the models
            maxi = 5;
            
            for modidx = maxi
                covfunc = {'covSum', {models{modidx}, 'covNoise'}}; % Set up covariance function - polynomial basis function of degree k with additive noise
                
                k = modidx ;                     % Degree of polynomial to fit (or if k == 5, model index of similarity function)
                if modidx < 5
                    logtheta = [-ones(1, k+1), -1]'; % Initialize the log hyperparameters
                    M = numel(logtheta);             % Nubmer of parameters
                else
                    logtheta = [-1, -1, -1]'; % Initialize the log hyperparameters
                    M = numel(logtheta);             % Nubmer of parameters
                end
                
                % Use gradient descent to find the optimal hyperparameters
                parms{i, modidx} = minimize(logtheta, 'gpr', -100, covfunc, xSTAR, ySTAR);
                
                [nlml, dnml, H] = gpr(parms{i,modidx}, covfunc, xSTAR, ySTAR);  % Compute the fit of the model - p(y|M) (marginal including hyperparameters), the partial derivatives and the Hessian matrix
                [fmu FS2] = gpr(parms{i,modidx}, covfunc, xSTAR, ySTAR, xPRED); % Compute the mean function and standard deviation for prediction space
                fS2 = FS2 - exp(2 * parms{i,modidx}(end));                      % Subtract off additive noise
                
                if isempty(strmatch('covPoly', models{modidx})) % if covSEiso
                    eparms = [exp(parms{i, modidx}(1)), exp(2 * parms{i, modidx}(2)), exp(2 *parms{i, modidx}(3))]';
                else
                    eparms  = [exp(-2 * parms{i, modidx}(1:(M - 1))); exp(2 * parms{i, modidx}(M))]'; % Exponentiated hyperprior means
                end

                K = feval(covfunc{:}, parms{i, modidx}, X);
                
                fullfits.means{i,modidx}  = fmu;
                fullfits.vars{i,modidx}   = fS2;
                fullfits.nParms(i,modidx) = numel(parms{i, modidx});
                fullfits.fits(i,modidx)   = -nlml; % Negative Log Marginal likelihood (higher values are better)
                
                S = []; 
                lnPrior = 0;
                if isempty(strmatch('covPoly', models{modidx})) % covSEiso

                    
                    % Length scale 
                    syms x m s positive
                    f = 1/(sqrt(2 * pi * s^2)) * exp(-1/(2 * s^2) * (x-m)^2); % Normal Distribution
                    h = -diff(log(f), x, 2); % Second derivative w.r.t. x
                    
                    pm = .55; ps = .15; % Prior parameters for length scale
                    lnPrior = lnPrior + subs(log(f), [x, m, s], [eparms(1), pm, ps]);
                    S(1, 1) = subs(h, [x, m, s], [eparms(1), pm, ps]);
                    
                    % Signal variance
                    f = 1/(sqrt(2 * pi * s^2)) * exp(-1/(2 * s^2) * (x-m)^2); % Normal Distribution
                    h = -diff(log(f), x, 2); % Second derivative w.r.t. x
                    
                    pm = 50; ps = 20; % Prior parameters for signal scale
                    lnPrior = lnPrior + subs(log(f), [x, m, s], [eparms(2), pm, ps]);
                    S(2, 2) = subs(h, [x, m, s], [eparms(2), pm, ps]);          
                        
                    % Parameters for similarity noise
                    pm = .08;
                    ps = .023;
                    parmIdx = 2; 
                else % covPoly
                    % Polynomial parameters
                    for parmIdx = 1:(M-1)
                        syms x m s positive  
                        f = 1/(sqrt(2 * pi * s^2)) * exp(-1/(2 * s^2) * (x-m)^2); % Normal Distribution
                        h = -diff(log(f), x, 2); % Second derivative w.r.t. x

                        pm = .05;  % Mean of polynomial parameter prior is 0
                        ps = .02; % Standard deviation of polynomial parameter prior is 10
                        lnPrior = lnPrior + subs(log(f), [x, m, s], [eparms(parmIdx), pm, ps]);
                        S(parmIdx, parmIdx) = subs(h, [x, m, s], [eparms(parmIdx), pm, ps]);
                    end
                    
                    % Parameters for polynomial noise
                    pm = .08;
                    ps = .023;
                end
                
                f = 1/(sqrt(2 * pi * s^2)) * exp(-1/(2 * s^2) * (x-m)^2); % Normal Distribution
                h = -diff(log(f), x, 2); % Second derivative w.r.t. x
                lnPrior = lnPrior + subs(log(f), [x, m, s], [eparms(parmIdx+1), pm, ps]);
                S(parmIdx+1, parmIdx+1) = subs(h, [x, m, s], [eparms(parmIdx+1), pm, ps]);

                A = -(H + S); 
                fullfits.lpD(i,modidx) = -nlml + lnPrior + .5 * M * log(2 * pi) - .5 * real(log(det(A))); % Laplace approximation to pD (Bishop p. 233)
                fullfits.lpDc(i,modidx) = -nlml + lnPrior + .5 * M * log(2 * pi) - .5 * log(det(A)); % Laplace approximation to pD (Bishop p. 233)
                % 8/5/2013 - Hack to avoid imaginary numbers
                
                % These are the normal BIC method divided by 2; the
                % ordinary method is -2 * lnml + M * log(numel(xSTAR))
                fullfits.bic(i,modidx) =     -nlml - .5 * M * log(numel(xSTAR)); % Larger numbers are better (Bishop2006)
                fullfits.bic2(i,modidx) = 2 * nlml +      M * log(numel(xSTAR));  % Smaller numbers are better (Wagenmakers2004)
                fullfits.draper(i,modidx) = -nlml     - .5 * M * log(numel(xSTAR)) + .5 * M * log(2 * pi); % Refer Chickering1996, the Draper measure is better than BIC
                fullfits.draper2(i,modidx) = 2 * nlml +      M * log(numel(xSTAR)) -      M * log(2 * pi); % Refer Chickering1996, the Draper measure is better than BIC
            end
            % Note that the bicWeights == bicWeights2 if the prior model
            % probabilities are equal (bicWeights treats bic as p(D|M) and computes the posterior p(M|D))
            fullfits.bicWeights(i,:)  = (exp(fullfits.bic(i, :)) .* exp(logmodelprior))/(exp(fullfits.bic(i,:)) * exp(logmodelprior)');
            fullfits.bicWeights2(i,:)  = exp(-.5 * fullfits.bic2(i,:) - min(fullfits.bic2(i,:)))./sum(exp(-.5 * fullfits.bic2(i,:) - min(fullfits.bic2(i,:))));
            
            fullfits.draperWeights(i,:)  = (exp(fullfits.draper(i, :)) .* exp(logmodelprior))/(exp(fullfits.draper(i,:)) * exp(logmodelprior)');
            fullfits.draperWeights2(i,:) = exp(-.5 * fullfits.draper2(i,:) - min(fullfits.draper2(i,:)))./sum(exp(-.5 * fullfits.draper2(i,:) - min(fullfits.draper2(i,:))));
            fullfits.pM(i,:)     = (exp(fullfits.lpD(i, :)) .* exp(logmodelprior))/(exp(fullfits.lpD(i,:)) * exp(logmodelprior)');
            
        end
        fullfits.X{i}     = X;
        fullfits.Y{i}     = Y;
        fullfits.xSTAR{i} = xSTAR;
        fullfits.ySTAR{i} = ySTAR;
    end
    fullfits.parms       = parms;
    toc
    save(fullfile(pwd, sprintf('gpAnalysis_Full_s%d.mat', subjects(sidx2))))
end