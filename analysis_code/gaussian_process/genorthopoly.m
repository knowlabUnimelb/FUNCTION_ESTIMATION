function pl = genorthopoly(x, n)
% Get all Chebeyshev coefficients from 0 to k
% This is my edit of orthogonalPoly - this is faster

kf = 1;

% This is a code downloaded from the website of MIT.
% http://ceta.mit.edu/comp_spec_func/

%       ==========================================================
%       Purpose: Compute orthogonal polynomials: Tn(x) or Un(x),
%                or Ln(x) or Hn(x), and their derivatives
%       Input :  KF --- Function code
%                       KF=1 for Chebyshev polynomial (First kind) Tn(x)
%                       KF=2 for Chebyshev polynomial (Second kind) Un(x)
%                n ---  Order of orthogonal polynomials
%                x ---  Argument of orthogonal polynomials
%       Output:  PL(n) --- Tn(x) or Un(x) or Ln(x) or Hn(x)
%                DPL(n)--- Tn'(x) or Un'(x) or Ln'(x) or Hn'(x)
%       =========================================================

% The only improvement in this program is it accepts vector arguments for x

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
lenx = length(x);

a = 2; b = 0; c = 1;
y0 = 1;
y1 = x;

pl = ones(lenx, n + 1);
if n > 0
    pl(:,2) = x; % T(1)

    if n > 1
        for  k = 3:(n + 1)
            yn = (a .* x + b) .* y1 - c * y0;
            pl(:,k) = yn;
            y0 = y1;
            y1 = yn;
        end
    end
end
if rowvec
    pl = pl';
end