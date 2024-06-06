function parms = dpgmm(data, varargin)
% function parms = dpgmm(data, num_its, parms)
% Dirichlet process Gaussian mixture model
% "rao-blackwellised" form, which does not store explicit means or covs
%
% If you supply parms as a variable, the count will start where the last iteration left off
%
% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt
%
% Edited by Daniel Little 2011
%% Compute optional variable
[T dimension] = size(data);

%some stats
debug = false;
allmean = mean(data,1);
allcov  = cov(data);

%% Assign optional variable
% a = T/50;
a = T/100;
optargs = {1000, 5, 100, [], a, .01, 6, allmean, allcov/10, 0, 0, [], [], ones(T, 1)};

% skip any new inputs if they are empty
newVals = cellfun(@(x) ~isempty(x), varargin);

% now put these defaults into the valuesToUse cell array, and overwrite the ones specified in varargin.
optargs(newVals) = varargin(newVals);
% or ... [optargs{1:numvarargs}] = varargin{:};

% Place optional args in memorable variable names
[nsamples, thinning, burnin, parms, alpha, kappa, nu, initmean, initcov, num_classes, counts, sums, cholSSE, classes] = optargs{:};


%% Set up parms structure)
alpha
if size(parms, 2) <= 1
    parms(1).alpha       = alpha;     % 1 / wishrnd(1,1); Prior probability of joining a new group
    parms(1).kappa       = kappa;     % A pseudo count on the mean (i.e., the prior mean is based on kappa pseudo-observations) [.1, T/1000]
    parms(1).nu          = nu;        % A pseudo-count on the covariance [6]
    parms(1).initmean    = initmean;  % Initial mean      [allmean]
    parms(1).initcov     = initcov;   % Intial covariance [allcov/10]
    parms(1).num_classes = num_classes;            % Initial number of classes
    parms(1).counts      = counts;            % Initial counts
    parms(1).sums        = sums;           % Initial sums
    parms(1).cholSSE     = cholSSE;           % Initial cholSSE
    parms(1).classes     = classes;    % Initial group assignment
    parms(1) = addNewClass(parms(1));                 %
    parms(1) = unhideObservations(parms(1), 1, data); %
    if debug, if ~checkParams (parms(1), data), disp('no check'); end, end
    continuedSampling = false;
else
    % Do nothing and use existing parms structure
    continuedSampling = true;
end

num_its = (nsamples * thinning) + burnin; 

%%
start_it = 1 + size(parms, 2); % Current count
for it = start_it:(start_it + num_its - 1)
    parms(it) = parms(it - 1); % Current parms equal previous parms
    
    % Disp iteration count, current number of clusters and current counts
    if mod(it - 1, 100) + 1 == 100
        disp(strcat(sprintf('%i [%i]: ', it, parms(it).num_classes), sprintf(' %i',parms(it).counts)));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GIBBS SAMPLING %%%%%%%%%%%%%%%%%%%%%%%
    t_order = randperm(T);
    for sctr = 1:T
        %t = sctr+1;%for debugging
        t = t_order(sctr); % Cycle observations in a random order
        
        old_class = parms(it).classes(t); % Retrieve old class assignment
        
        parms(it) = hideObservations(parms(it), old_class, data(t,:)); % Update parms for clusters after removing the current item
        parms(it) = handleRemovedClasses(parms(it));                   % Deal with clusters that disappear
        if debug, if ~checkParams(parms(it),data,t), disp('no check at hide'); end, end
        
        % these are the probabilities that we will sample from
        % note we add one to include the chance of adding a new class
        log_p_obs = -inf * ones(parms(it).num_classes + 1, 1); % Initialize probabilities for joining each group
        
        p_prior = [];
        parms(it) = addNewClass(parms(it));  % add a new cluster to the dpgmm, it will be removed if nothing is put in it
        if debug, if ~checkParams(parms(it),data,t), disp('no check at add class'); end, end
        
        kappabar = parms(it).counts + parms(it).kappa; % Posterior kappa (count for each group)
        nubar    = parms(it).counts + parms(it).nu;    % Posterior nu
        factor   = (kappabar + 1) ./ (kappabar .* (nubar - dimension - 1));
        p_prior  = parms(it).counts + parms(it).alpha * (parms(it).counts == 0);
        
        if dimension == 1
            log_p_obs = mvnormpdfln(data(t,:)', ...    % Current observation
                parms(it).sums' ./ kappabar,...        % Mean for each cluster
                (sqrt(factor)' .* squeeze(parms(it).cholSSE))')'; % Std for each cluster
        else
            for i = 1:parms(it).num_classes
                %            if (parms(it).counts(i) == 0), p_prior(i) = parms(it).alpha;
                %            else p_prior(i) = parms(it).counts(i); end
                try
                    % integrating over the parameters of a normal-inverse-Wishart yields student-t.
                    % This can be approximated by a "moment-matched" Gaussian, see Erik Sudderth (2006) Thesis, p. 47
                    % kappabar = parms(it).counts(i) + parms(it).kappa;
                    % nubar = parms(it).counts(i) + parms(it).nu;
                    % factor = (kappabar + 1) / (kappabar * (nubar - dimension - 1));
                    log_p_obs(i) = mvnormpdfln(data(t,:)', ...        % Current observation
                        parms(it).sums(i,:)' / kappabar(i),...        % Mean
                        sqrt(factor(i)) * parms(it).cholSSE(:,:,i));
                catch
                    disp('mvnpdf throws error');
                end
            end
        end
        
        %lightspeed sample normalizes automatically
        classprobs = p_prior' .* exp(log_p_obs - max(log_p_obs));
        try
            new_class = sample(classprobs);
            %             if (parms(it).counts(new_class) == 0)
            %                disp('adding a guy');
            %             end
            parms(it).classes(t) = new_class;
        catch
            disp('could not sample');
        end
        
        if debug, if ~checkParams(parms(it),data,t), disp('no check at sample'); end, end
        
        parms(it) = unhideObservations(parms(it), new_class, data(t,:));
        if debug, if ~checkParams(parms(it),data), disp('no check at hide'); end, end
        
    end
    
    %%%%%%%%%%%%%%%%%%%% PARAMETER UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %alpha is the "pseudo-count" for new classes.  it is estimated using ARS
    k = parms(it).num_classes;
    n = T;
    
    %can show that derivative is guaranteed to be positive / negative at
    %these points
    deriv_up = 2 / (n - k + 3/2);
    deriv_down = k * n / (n - k + 1);
    
    % this is the version with a conjugate inverse gamma prior on alpha, as
    % in Rasmussen 2000
    parms(it).alpha = ars(@logalphapdf, {k, n}, 1, [deriv_up deriv_down], [deriv_up inf]);
    
    %this is the version with a totally non-informative prior
%     parms(it).alpha = ars(@logalphapdfNI, {k, n}, 1, [deriv_up deriv_down], [deriv_up inf]);
end

%% Sampling can not be continued using the code below
if ~continuedSampling
    parmIdx = 1:numel(parms);
    parmIdx(1:burnin) = [];
    pi = parmIdx(1:thinning:(nsamples * thinning));
    for i = 1:numel(pi)
        outparms(i) = parms(pi(i));
    end
    parms = outparms;
end

%%
%checks a set of parameters to see if they are self-consistent.
%for debugging
function [total c_basic c_count c_sum] = checkParams(parms,data,exclude)
if exist('exclude','var')
    c_basic = min(parms.classes([1:exclude-1 exclude+1:end]) > 0);
else
    c_basic = min(parms.classes > 0);
end
c_count = 1;
c_sum = 1;
for i = 1:parms.num_classes
    statedata = data(find(parms.classes == i),:);
    err_amount = parms.sums(i,:) - sum(statedata) - parms.kappa * parms.initmean;
    statecount = size(statedata,1);
    if exist('exclude','var')
        if i == parms.classes(exclude)
            err_amount = err_amount - data(exclude,:);
            statecount = statecount - 1;
        end
    end
    if (statecount ~= parms.counts(i)), c_count = 0; end
    if (sum(err_amount) > .01), c_sum = 0; end
end
total = c_basic * c_count * c_sum;