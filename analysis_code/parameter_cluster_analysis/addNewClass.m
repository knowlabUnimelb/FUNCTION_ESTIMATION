% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

function parms = addNewClass(parms)
%function parms = addNewClass(parms)
%adds a new, empty class to the dpmm

newclassidx        = parms.num_classes + 1;
parms.num_classes = newclassidx;
parms.counts(newclassidx) = 0;
parms.sums(newclassidx,:) = parms.kappa * parms.initmean;
parms.cholSSE(:,:,newclassidx) = chol(parms.nu * parms.initcov);
%parms.SSE(:,:,newclassidx) = parms.nu * parms.initcov;
