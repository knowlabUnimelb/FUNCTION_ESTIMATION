% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt
%
% Updates the Normal-Inverse-Wishart Parameters

function [parms oldparms] = hideObservations(parms, new_class, data)%
% updates parms when a given set of observations are hidden.  
% this is necessary for Gibbs sampling
oldparms = parms;

% update the suff stats see Sudderth's (1999) thesis (page 46)

% Get current cluster count for selected data point
old_count = parms.counts(new_class) + parms.kappa; % Add kappa so that the count is the likelihood count + the prior count
old_sum   = parms.sums(new_class,:); % Get previous cluster sum

% Update the posterior hyperparameters by cacheing the observations’ sum (
%   eq. (2.62)), and the Cholesky decomposition
%   of the sum of observation outer products (eq. (2.63)). Cholesky decompositions 
%   are numerically robust, can be recursively updated as observations are added or removed, and
%   allow fast likelihood evaluation through the solution of triangulated linear systems.
parms.cholSSE(:,:,new_class) = cholupdate(parms.cholSSE(:,:,new_class), old_sum' / sqrt(old_count));
% parms.SSE(:,:,new_class) = parms.SSE(:,:,new_class) + old_sum' * old_sum / old_count;

% must add iteratively because of the cholesky update function
for i = 1:size(data,1)
    
    % Subtract current data point from the count of its old cluster
    parms.counts(new_class) = parms.counts(new_class) - 1;
    % Subtract data from sum
    parms.sums(new_class,:) = parms.sums(new_class,:) - data(i,:);
   
    try
        parms.cholSSE(:,:,new_class) = cholupdate(parms.cholSSE(:,:,new_class), data(i,:)','-');
    catch
        disp('Error in hideObservations at 33. CHOLUPDATE failed.');
    end
    %   parms.SSE(:,:,new_class) = parms.SSE(:,:,new_class) - data(i,:)' * data(i,:);
end


new_count = parms.counts(new_class) + parms.kappa;
parms.cholSSE(:,:,new_class) = cholupdate(parms.cholSSE(:,:,new_class), parms.sums(new_class,:)' / sqrt(new_count), '-');
%parms.SSE(:,:,new_class) = parms.SSE(:,:,new_class) - new_sum' * new_sum / new_count;