function [mrcnt, rd] = plotDrawResponse(drawdata, nr, dataSets, firstline, i, ploton)
% global ploton
% drawdata = full dataset
% nr = number of responses
% dataSets = data set IDs
% firstline = firstline of data set
% i = set count

% Cycle through the response data for the normal dataset
mrcnt = [];
rd = [];
for ri = 1:numel(nr)
    if ri == 2 % Tag trialnumber for later
        mrcnt = [mrcnt; dataSets(i) firstline(1,2) firstline(1,12)]; % data set ID, trial number, 1 = normal/2 = reversed
    end
    responsedata = drawdata(drawdata(:,2) == firstline(1,2) & drawdata(:,3) == nr(ri), :); % Get only response ri
    
    if ploton
%         if responsedata(1,12) == 1
            plot(responsedata(responsedata(:,6) ~= 0,6), responsedata(responsedata(:,6) ~= 0,7), 'LineWidth', 1.5) % Plot that response in blue
%             plot(1 * responsedata(responsedata(:,6) ~= 0,6), 1 * responsedata(responsedata(:,6) ~= 0,7), 'LineWidth', 1.5); hold on % Plot that response in blue
%         elseif responsedata(1,12) == 2
%             plot(-1 * responsedata(responsedata(:,6) ~= 0,6),  -1 * responsedata(responsedata(:,6) ~= 0,7), 'g', 'LineWidth', 1.5) % Plot that response reversed in green
% %             plot(-1 * responsedata(responsedata(:,6) ~= 0,6),  1 * responsedata(responsedata(:,6) ~= 0,7), 'g', 'LineWidth', 1.5); hold on % Plot that response reversed in green
%         end
    end
    rd = [rd; responsedata];
end