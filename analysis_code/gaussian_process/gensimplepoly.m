function pl = gensimplepoly(x, n, bias)
% Get all Chebeyshev coefficients from 0 to k
% This is my edit of orthogonalPoly - this is faster

if nargin == 2
    bias = 1;
end

% make sure that x is a row or column vector and not a matrix.
[r, c] = size(x);
if r == 1
    x = x';
    rowvec = 1;
elseif c == 1
    rowvec = 0;
else
    error('x must be a vector, and cannot be a matrix');
end

if bias == 1
for i = 0:n
    pl(:,i+1) = x .^ i;
end
else
    for i = 1:n
    pl(:,i+1) = x .^ i;
end
end

if rowvec
    pl = pl';
end