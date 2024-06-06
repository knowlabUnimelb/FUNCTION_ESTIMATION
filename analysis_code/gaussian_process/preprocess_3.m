% Script to plot FUNCTION ESTIMATION data from same subjects to
% examine reliability of drawing method
%
% This method uses the final_s*_data.mat file that contains the data
% corrected for multiple responses and sorted
%
% The output is a mat file containing finalData matrices of sorted data,
% containing only one response with the noise data adjusted to the
% appropriate size and the draw data replicated so that there are only
% nan's where there are no observed or response data - all other trial
% information is replicated to make indexing easier
%
% The initial function was written for subjects who participate in both
% the mouse version and the drawing version. This version is for the
% drawing data only (16/8/10)
function preprocess_3(subject, ploton, saveon)
clc
close all


%% Load drawing data
if nargin == 0
    subject = 301;
    ploton = true;
elseif nargin == 1
    ploton = true;
end

datafolder = 'data';

pploton = ploton;
if exist(fullfile(pwd, datafolder, sprintf('multiple_resp_corrected_s%d_data.mat', subject)), 'file') == 2
    load(fullfile(pwd, datafolder, sprintf('multiple_resp_corrected_s%d_data.mat', subject)));
else
    load(fullfile(pwd, datafolder, sprintf('preprocess1_out_s%d_data.mat', subject)));
end
ploton = pploton;

dataSets = 1:26;
noutliers = [0 1 3];
finalDrawData = [];

%%
% sets = {1:18; 19:36; 37:54; 55:72; 73:90; 91:96};
sets = {1:26};
for j = 1:size(sets, 1)
    if ploton; figure('WindowStyle', 'docked'); end
    pcnt = 1;
    for i = sets{j}
        % Plot Draw Data
        if ploton; subplot(5,6,pcnt); end
%         if subject < 901
            r1 = sortedDrawData(sortedDrawData(:,2) == dataSets(i),:);% & sortedDrawData(:,12) == 1, :);
%             r2 = sortedDrawData(sortedDrawData(:,13) == dataSets(i) & sortedDrawData(:,12) == 2, :);
%         else
%             r1 = [];
%             r2 = [];
%         end
        
        % Check if either are empty
        is1 = ~isempty(r1);
%         is2 = ~isempty(r2);
%         [i is1 is2]
        
        if is1 %&& is2
            d1 = sortedDrawData(sortedDrawData(:,2) == r1(1,2) & sortedDrawData(:,3) == 1, :);
%             d2 = sortedDrawData(sortedDrawData(:,2) == r2(1,2) & sortedDrawData(:,3) == 1, :);
            
            if ploton
%                 scatter(d1(:,4), -1 *  d1(:,5), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
                scatter(d1(:,4), d1(:,5), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
%                 plot(d1(d1(:,6) ~= 0,6), -1 * d1(d1(:,6) ~= 0,7), 'LineWidth', 1.5)
                plot(d1(d1(:,6) ~= 0,6), d1(d1(:,6) ~= 0,7), 'LineWidth', 1.5)
%                 plot(-1 * d2(d2(:,6) ~= 0,6),  -1 * d2(d2(:,6) ~= 0,7), 'g', 'LineWidth', 1.5)
            end
            
            finalDrawData = [finalDrawData; [d1(:,1:7), repmat(d1(1,8:13), size(d1(:,8:13), 1), 1)]];%;...
%                 [d2(:,1:7), repmat(d2(1,8:13), size(d2(:,8:13), 1), 1)]] ;
%         elseif is1 && ~is2
%             d1 = sortedDrawData(sortedDrawData(:,2) == r1(1,2) & sortedDrawData(:,3) == 1, :);
%             
%             if ploton
%                 scatter(d1(:,4),  -1 * d1(:,5), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
%                 plot(d1(d1(:,6) ~= 0,6),  -1 * d1(d1(:,6) ~= 0,7), 'LineWidth', 1.5)
%             end
%             
%             finalDrawData = [finalDrawData; [d1(:,1:7), repmat(d1(1,8:13), size(d1(:,8:13), 1), 1)]];
%         elseif ~is1 && is2
%             d2 = sortedDrawData(sortedDrawData(:,2) == r2(1,2) & sortedDrawData(:,3) == 1, :);
%             
%             if ploton
%                 scatter(-1 * d2(:,4),  -1 * d2(:,5), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
%                 plot(-1 * d2(d2(:,6) ~= 0,6),  -1 * d2(d2(:,6) ~= 0,7), 'g', 'LineWidth', 1.5)
%             end
%             
%             finalDrawData = [finalDrawData; [d2(:,1:7), repmat(d2(1,8:13), size(d2(:,8:13), 1), 1)]];
        end
        if ploton
            axis([-1 1 -1 1]); box on
            set(gca,'XTick', [], 'YTick', [])
            title(num2str(i))
            set(gca,'FontSize', 8)
        end
        pcnt = pcnt + 1;
    end
end

if saveon 
finalDrawData = finalDrawData([[~isnan(finalDrawData(:,4)) & ~isnan(finalDrawData(:,5))] | [~isnan(finalDrawData(:,6)) & ~isnan(finalDrawData(:,7))]], :);
 filestr = fullfile(pwd, datafolder, sprintf('final_%d_data', subject));
  eval(['save ''' filestr ''''])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%