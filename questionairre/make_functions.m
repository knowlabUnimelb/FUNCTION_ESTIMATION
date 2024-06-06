% For RUN_FUNCTION_NOISE, functions should be generated upfront and the same
% functions should be used for all experiments
%
% In this version of the function, there is a gap at the outlier position.
% The outlier positions are chosen arbitrarily. The full design is 16
% functions x 1 noise level x 2 number of data points x 4 outliers (none,
% near, far, three)
%
% NOTE: This looks good so far but I want to present exactly the same data
% set for each outlierCode position. So I need to make the functions and
% then duplicate them 4 times.
%
clear all
clc

%% Experimental Functions
% Load function information from loadFunc.txt
%   loadFunc.txt --> [polynomial_coefficients, noise(.25, .50), n(5, 15, 50), outlierLocation(0,1,2), includeIdx(1 if to be used, 0 if non-shown complement)]
%   there are positive and negative instances of each function degree, each
%   instantiation of type x noise x n x outlier is only shown once in
%   either a positive or negative orientation
fid = fopen('loadFunc.txt');
columnNames = textscan(fid, '%s %s %s %s %s %s %s %s', 1);
functions = cell2mat(textscan(fid, '%f %f %f %f %f %f %f %f'));
fclose(fid);

% Because I'm removing a 3rd of the points, I need to produce an extra 3rd
% of the points
functions(functions(:,7) == 5,7) = 10; 
functions(functions(:,7) == 25,7) = 40;
functions(:,8) = [];

functions = [functions; functions];
% sflaglocs = [-6 -2; -2  2; 2 6]; 
             
cnt = 1;
for i = 1:size(functions, 1)
    x{cnt} = linspace(-6, 6, functions(i,7)); % X values = linspace(xmin, xmax, n)
%     showflag{cnt} = ~((x{cnt} >= sflaglocs(functions(i,8),1)) & (x{cnt} <= sflaglocs(functions(i,8),2)));  % Select items to show based on gap location
    xhat{cnt} = linspace(-6, 6, 100); % X values for true x functions
    
    y{cnt} = gety(x{cnt}, functions(i,1:5), functions(i,6)); % Get yvalues from polynomial coefficients and x values
    yhat{cnt} = gety(xhat{cnt}, functions(i,1:5), 0); % Get yvalues from polynomial coefficients and x values

    x{cnt} = x{cnt} + normrnd(0, .25, size(x{cnt}, 1), size(x{cnt}, 2)); % Add horizontal jitter
%   UNCOMMENT TO CHECK FIGURE ACTION
%     figure
%     scatter(x{cnt}, y{cnt}); hold on
%     plot(xhat{cnt}, yhat{cnt}, '.k')
%     pause
%     close(gcf)
    
    cnt = cnt + 1;
    
end 
nuniquefun = cnt - 1;


trial_sequence = functions; %repmat(functions, 3, 1);
% x = repmat(x, 1, 3);
% y = repmat(y, 1, 3);
% xhat = repmat(xhat, 1, 3);
% yhat = repmat(yhat, 1, 3);
% showflag = repmat(showflag, 1, 3);

% trial_sequence(:,9) = reshape(repmat([0 1 3], nuniquefun, 1), nuniquefun * 3, 1);
% to = randperm(size(trial_sequence,1)); % Randomize trial order
% trial_sequence = trial_sequence(to,:);
% 
% x = x(to); y = y(to); xhat = xhat(to); yhat = yhat(to); showflag = showflag(to);

% mlocs = (sflaglocs(:,2) + sflaglocs(:,1))/2;
% cnt = 1;
% for i = 1:size(trial_sequence, 1)
%        % Generate outliers
%     if trial_sequence(i,9) == 1
%         % Query x and y points and add a point at the maximum  distance
%         qx = normrnd(mlocs(trial_sequence(i,8)), .5);
%         qy = linspace(min(y{cnt}), max(y{cnt}), 100);
%         qxy = allcomb(qx, qy);
%         
%         [minxhatval, minxhatidx] = min(abs(xhat{cnt} - qx));
%         dists = [];
%         dpos = repmat([xhat{cnt}(minxhatidx) yhat{cnt}(minxhatidx)], size(qxy,1), 1);
%         dists = sqrt((dpos(:,2) - qxy(:,2)).^2 + (dpos(:,1) - qxy(:,1)).^2);
%         [yy, ii] = max(dists);
%         
%         outlier{cnt} = qxy(ii,:);
%         distance{cnt} = yy;
%     elseif trial_sequence(i,9) == 3
%         % Query x and y points and add 3 points at the maximum  distance
%         qx = normrnd(mlocs(trial_sequence(i,8)), .5);
%         qy = linspace(min(y{cnt}), max(y{cnt}), 100);
%         qxy = allcomb(qx, qy);
% 
%         [minxhatval, minxhatidx] = min(abs(xhat{cnt} - qx));
%         dists = [];
%         dpos = repmat([xhat{cnt}(minxhatidx) yhat{cnt}(minxhatidx)], size(qxy,1), 1);
%         dists = sqrt((dpos(:,2) - qxy(:,2)).^2 + (dpos(:,1) - qxy(:,1)).^2);
%         [yy, ii] = max(dists);
%         
%         outlier{cnt}(:,1) = normrnd(qxy(ii,1), std(x{cnt})/10, 3, 1);
%         outlier{cnt}(:,2) = normrnd(qxy(ii,2), std(y{cnt})/10, 3, 1);
%         distance{cnt} = yy;
%     else
%         outlier{cnt} = [NaN NaN];
%     end 
%     
%         figure('WindowStyle', 'docked')
%     scatter(x{cnt}(showflag{cnt}), y{cnt}(showflag{cnt})); hold on
%     a = scatter(outlier{cnt}(:,1), outlier{cnt}(:,2));
%     set(a, 'MarkerFaceColor', [1 0 0]);
%     plot(xhat{cnt}, yhat{cnt}, '.k')
%     title(num2str(trial_sequence(i,9)))
%     pause
%     close(gcf)
% 
%     %   UNCOMMENT TO CHECK FIGURE ACTION
%     
%     cnt = cnt + 1;
% end