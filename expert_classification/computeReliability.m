clear all
clc

load handRatings
% x is rater 1
% y is rater 2
[r, p] = corrcoef(x(:), y(:), 'rows', 'complete');

% Cohen's Kappa
pa = nansum(x(:) == y(:))./numel(x) ;
pe = 1/max(x(:));
K = (pa - pe)./(1 - pe);
% 0–0.20 as slight, 0.21–0.40 as fair, 0.41–0.60 as moderate, 0.61–0.80 as substantial, and 0.81–1 as almost perfect agreement. 