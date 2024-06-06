% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

function parms = handleRemovedClasses(parms)

% This looks like you can set a certain number of clusters as fixed
if (~isfield(parms,'num_fixed') || parms.num_fixed == inf) 
    idxs = find(parms(end).counts == 0);
    for ctr = idxs
        
        % reduce all state numbers that are greater than ctr
        % If a cluster is removed, then shift all of the other cluster
        % labels down by 1, if the cluster is at the end, do nothing
        parms.classes = parms.classes - (parms.classes >= ctr); 
        
        % Get new cluster labels
        idxs2 = [1:(ctr-1) (ctr+1):parms.num_classes];
        
        
        parms.counts = parms.counts(idxs2); % Adjust counts to correct for removed cluster
        parms.sums = parms.sums(idxs2,:);   % Adjust sum to correct for removed cluster
        parms.cholSSE = parms.cholSSE(:,:,idxs2);   % Adjust cholSSE to correct for removed cluster
        %        parms.SSE = parms.SSE(:,:,idxs);
        parms.num_classes = parms.num_classes - 1;  % Adjust num_classes to correct for removed cluster
    end
end