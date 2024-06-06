% Script to take a sorted_s*_data.mat file, display trials with multiple
% responses, select or combine multiple responses and save
%
% The output is a matfile containing data with the multiple responses
% corrected so that there is only one response per trial
%
% This script was initially written for an analysis looking at subjects who
% completed both the drawing and the mouse version of the experiment. As of
% 16/8/2010, I've edited the program to only examing the drawing version of
% the experiment
function preprocess_2_correct_multiple_responses(subject)
close all

datafolder = 'data';

if nargin == 0
    subject = 301; 
end

if exist(fullfile(pwd, datafolder, sprintf('multiple_resp_corrected_s%d_data.mat', subject)), 'file') == 2
    load(fullfile(pwd, datafolder, sprintf('multiple_resp_corrected_s%d_data.mat', subject)));
elseif exist(fullfile(pwd, datafolder, sprintf('preprocess1_out_s%d_data.mat', subject)), 'file') == 2
    load(fullfile(pwd, datafolder, sprintf('preprocess1_out_s%d_data.mat', subject)))
    currtrial = 1;  
else
    error('File not found')
end
% currtrial = 83;
% colors = jet; set(gcf, 'WindowStyle', 'docked')
disp(['Total trials to adjust = ' num2str(size(mrcnt,1))])
for idx = currtrial:size(mrcnt, 1) % Cycle through problems with multiple responses
    close all
    figure('WindowStyle', 'docked')
    
    % Get data from currently selected trial
    % mrcnt = [dataset, trial number, normal/reversed]
    cd = sortedDrawData(sortedDrawData(:,2) == mrcnt(idx, 2), :);
    cdStartIdx = find(all(sortedDrawData == repmat(cd(1,:), size(sortedDrawData, 1), 1), 2));
    
    befdata = sortedDrawData(1:(cdStartIdx-1),:);
    aftdata = sortedDrawData((cdStartIdx + size(cd,1)):end,:);
    
    % Plot current problem
    colors = jet(max(cd(:,3)));
    h = ones(max(cd(:,3)), 1);
    for ridx = 1:max(cd(:,3)) % cycle through each response
        rd = cd(cd(:,3) == ridx,:);
        
        coloridx = colors(round(linspace(1, size(colors, 1), max(cd(:,3)))), :);
        
        scatter(rd(:,4), -1 *  rd(:,5), 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
        %         if mrcnt(idx, 3) == 1
        h(ridx) = plot(rd(rd(:,6) ~= 0,6), -1 * rd(rd(:,6) ~= 0,7), 'Color', coloridx(ridx,:), 'LineWidth', 1.5); % Plot that response in blue
        %         else
        %             h(ridx) = plot(-1 * rd(rd(:,6) ~= 0,6),  -1 * rd(rd(:,6) ~= 0,7), 'g', 'LineWidth', 1.5); % Plot that response reversed in green
        %         end
    end
    legend(h, num2str((1:ridx)'), 'Location', 'BestOutside')
    axis([-1 1 -.8 .8]); box on
    
    %% Select which responses to keep or combine
    selected = input('Which functions do you want to keep?');
    
    if numel(selected) == 2 % If more than one is selected
        overlap  = input(sprintf(['Which functions overlap: [' repmat('%d ', 1, numel(selected)) ': 0 = none]?'], selected));
        if overlap == 0
            combine = 2; % No overlap - interpolate
        else
            combine  = 1;
        end
    elseif numel(selected) > 2
        combine = 2;
    else % If only one is selected
        overlap = [];
        combine = 1;
    end
    
    % Combine or select responses
    if combine == 1 % Just replace points at the overlapping points from the chosen function
        if numel(selected) > 1
            pick = [];
            while isempty(pick)
                pick = input('Which function do you want to keep?');
                if numel(pick) > 1
                    pick = [];
                    disp('Pick one function to keep')
                end
            end
        end
        
        sd = cd(ismember(cd(:,3), selected),:); % selected data
        if numel(selected) == 1 % If only one is selected
            outsd = cd(cd(:,3) == 1,:); nobs = size(outsd,1);
            if size(sd, 1) > size(outsd, 1)
                sizediff = size(sd, 1) - size(outsd, 1);
                outsd = [outsd; repmat(outsd(end,:), sizediff, 1)];
                 outsd((nobs + 1):end, 4:5) = nan;
                outsd(:,6:7) = sd(:,6:7);
            elseif size(sd, 1) < size(outsd, 1)
                sizediff = size(outsd, 1) - size(sd, 1);
                ssd = sd(:,6:7);
                ssd = [ssd; nan(sizediff, 2)];
                outsd(:,6:7) = ssd;
            else
                outsd(:,6:7) = sd(:,6:7);
            end
            sortedDrawData = [befdata; outsd; aftdata]; % Just keep that response
            out = cd(ismember(cd(:,3), selected),6:7);
        elseif numel(selected) == 2 % Find the overlapping section
            for sidx = 1:numel(selected)
                sresp{sidx} = sd(sd(:,3) == selected(sidx),6:7);
            end
            minpoint = min(sresp{pick == selected}(:,1));
            maxpoint = max(sresp{pick == selected}(:,1));
            
            if minpoint < min(sresp{pick ~= selected}(:,1))
                out = [];
                out = [out; sresp{selected == pick}];
                out = [out; sortrows([sresp{selected ~= pick}(sresp{selected ~= pick}(:,1) < minpoint, :);...
                    sresp{selected ~= pick}(sresp{selected ~= pick}(:,1) > maxpoint, :)], 1)];
            else
                out = [];
                out = [out; sortrows([sresp{selected ~= pick}(sresp{selected ~= pick}(:,1) < minpoint, :);...
                    sresp{selected ~= pick}(sresp{selected ~= pick}(:,1) > maxpoint, :)], 1)];
                out = [out; sresp{selected == pick}];
            end
            
            outsd = cd(cd(:,3) == 1,:); nobs = size(outsd,1);
            if size(out, 1) > size(outsd, 1)
                sizediff = size(out, 1) - size(outsd, 1);
                outsd = [outsd; repmat(outsd(end,:), sizediff, 1)];
                outsd((nobs + 1):end, 4:5) = nan;
                outsd(:,6:7) = out;
            elseif size(out, 1) < size(outsd, 1)
                sizediff = size(outsd, 1) - size(sd, 1);
                out = [out; nan(sizediff, 1)];
                outsd(:,6:7) = out;
            else
                outsd(:,6:7) = out;
            end
            sortedDrawData = [befdata; outsd; aftdata]; % Just keep that response
        end
    elseif combine == 2
        ok2sort = input('Which of these are ok to sort?');
        out = [];
        for ssi = 1:numel(selected)
            if ismember(ssi, ok2sort)
                out = [out; sortrows(cd(ismember(cd(:,3), selected(ssi)),6:7), 1)];
            else
                out = [out; (cd(ismember(cd(:,3), selected(ssi)),6:7))];
            end
        end
        outsd = cd(cd(:,3) == 1,:); nobs = size(outsd,1);
        if size(out, 1) > size(outsd, 1)
            sizediff = size(out, 1) - size(outsd, 1);
            outsd = [outsd; repmat(outsd(end,:), sizediff, 1)];
            outsd((nobs + 1):end, 4:5) = nan;
            outsd(:,6:7) = out;
        elseif size(out, 1) < size(outsd, 1)
            sizediff = size(outsd, 1) - size(out, 1);
            out = [out; nan(sizediff, 2)];
            outsd(:,6:7) = out;
        else
            outsd(:,6:7) = out;
        end
        sortedDrawData = [befdata; outsd; aftdata]; % Just keep that response
    end
    
    % Plot new data
    plot(out(:,1), -1 * out(:,2), 'Color', 'k', 'LineWidth', 2); % Plot that response in black
    axis([-1 1 -.8 .8]); box on
    pause
    
    % Save data file
    currtrial = currtrial + 1;
    filestr = fullfile(pwd, datafolder, sprintf('multiple_resp_corrected_s%d_data', subject));
    eval(['save ''' filestr ''''])
end
