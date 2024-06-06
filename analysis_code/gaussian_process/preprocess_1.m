% Script to plot FUNCTION ESTIMATION data from same subjects to
% examine reliability of drawing method
%
% The output is a matfile of sorted data (by dataset and normal/reversed)
% which can be passed to correct_multiple_responses.m
%
% This script was initially written for subjects who completed both the
% drawing and mouse version. As of 16/8/2010, this version has been modified
% to accomodate only the drawing version. This script will create a mat
% file which identifies the trials that have multiple responses. This mat
% file can then be passed to correct_multiple_responses.

function preprocess_1(subject, ploton, saveon)
close all

%% Load drawing data

datafolder = 'data';

if exist(fullfile(pwd, datafolder, sprintf('s%03d_socialFunctionData.dat', subject)), 'file') == 2
    drawdata = dlmread(fullfile(pwd, datafolder, sprintf('s%03d_socialFunctionData.dat', subject)));
else
    error(sprintf('Data file missing for subject %d', subject))
end

ntrials = max(drawdata(:,2));

dataSets = 1:26;

%%
sets = {1:26}; % Plot 13 separate plots per figure panel
mrcnt = []; % Initialize multiple response count
sortedDrawData  = []; % Initialize

for j = 1:size(sets, 1)
    if ploton; figure('WindowStyle', 'docked'); end % For each set of 18, open a new figure
    
    pcnt = 1; % Initialize subplot count
    for i = sets{j} % Cycle through current set
        if ploton; subplot(5, 6, pcnt); end

        %% GENDRAW
        % Find normal and reversed responses to the same data set (only the
        % first line of each entry will have all of the trial information)
        r1 = drawdata(drawdata(:,2) == dataSets(i), :); % column 12 = 1 if normal
%         r2 = drawdata(drawdata(:,13) == dataSets(i) & drawdata(:,12) == 2, :); % column 12 = 2 if reversed
        
        % Check if empty (indicating that the response was missed)
        is1 = ~isempty(r1);

        % Get number of actual responses
        if is1 % If there is a response
            % Now get all of the data for each problem
            d1 = drawdata(drawdata(:,2) == r1(1,2) & drawdata(:,3) == 1, :); % Normal
            
            nr1 =  unique(drawdata(drawdata(:,2) == r1(1,2), 3)); % column 3 indexes the response count
            
            % Plot the observed data
            if ploton
                scatter(d1(:,4), d1(:,5), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
            end
            
            % Plot the response from the normal data set and return the
            % multiple response count
            [multiResp, rd1] = plotDrawResponse(drawdata, nr1, dataSets, r1(1,:), i, ploton);
            mrcnt = [mrcnt; multiResp];
            
%             % Plot the response from the reversed data set and return the
%             % multiple response count
%             [multiResp, rd2] = plotDrawResponse(drawdata, nr2, dataSets, r2, i, ploton);
%             mrcnt = [mrcnt; multiResp];
            
            % Keep the sorted data in order by data set id
            sortedDrawData = [sortedDrawData; rd1]; %; rd2];
           
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
supertitle(sprintf('Subject %d', subject))
if saveon
filestr = fullfile(pwd, datafolder, sprintf('preprocess1_out_s%d_data', subject));
eval(['save ''' filestr ''''])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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