clear all
clc
close all

nProblems = 24; 

load zParms

Z = reshape(z, numel(z)/nProblems, nProblems);
NV = reshape(znv, numel(znv)/nProblems, nProblems);

currentdirectory = pwd;
cd ..
parentdirectory = pwd;
cd(currentdirectory)

datafiles = dir(fullfile(parentdirectory, 'Fits', '*.mat'));

lengthScale = nan(size(datafiles, 1), nProblems);
signalParm = nan(size(datafiles, 1), nProblems);
noiseParm = nan(size(datafiles, 1), nProblems);
for ii = 1:size(datafiles, 1)
    load(fullfile(parentdirectory, 'Fits', datafiles(ii).name));
%     [a,b] = strtok(datafiles(ii).name, '.');
%     subIdx(ii,1) = str2double(a(18:end));
    subjectFits(ii) = fullfits;
%     for jj = 1:nProblems
%         if ~isempty(fullfits.parms{jj,5})
%             for kk = 1:3
%                 covParms(ii,jj,kk) = fullfits.parms{jj,5}(kk);
%             end
%         end
%     end
end

%% Cluster using dpgmm
markers = {' xb', ' og', ' sr', ' ^c', ' +m', '*k'};
lineColors = {'b', 'g', 'r', 'c', 'm', 'k'};
minCutoffCount = 10; 
for i = 1%nProblems
    clear handle legNames
    figure('WindowStyle', 'docked'); 
    clusterData{i} = [Z(:,i), NV(:,i)]; 
    clusterOutput{i} = dpgmm(clusterData{i});
    
    clusters{i} = mode([clusterOutput{i}(:).classes], 2);
    
    counts{i} = aggregate([clusters{i}, ones(size(clusters{i},1),1)], 1, 2, @count);
    discardedcounts{i} = counts{i}(counts{i}(:,2) < minCutoffCount, :);
    counts{i}(counts{i}(:,2) < minCutoffCount, :) = [];
    
    
    for j = 1:size(counts{i}, 1)
         subplot(1,2,1)
         handle(j) = plot(clusterData{i}(clusters{i} == counts{i}(j,1),  1),... 
                           clusterData{i}(clusters{i} == counts{i}(j,1),  2),... 
                           markers{j});
         xlabel('log Length Scale + Signal Var', 'FontSize', 12)
         ylabel('log Noise Variance', 'FontSize', 12)
%          zlabel('log Noise Variance', 'FontSize', 12)
         hold on
         box on
         legNames{j} = sprintf('Cluster %d (N = %d)', j, counts{i}(j,2));
         
         subplot(1,2,2)

         
         clusterSubs = find(clusters{i} == j);
         for k = 1:numel(clusterSubs)
            plot(subjectFits(clusterSubs(k)).xSTAR{i},...
                 subjectFits(clusterSubs(k)).ySTAR{i},...
                 'Color', lineColors{j}, 'Marker', 'none')
             hold on
         end
         plot(subjectFits(1).X{i}, subjectFits(1).Y{i}, ' ok', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
         xlabel('X', 'FontSize', 12)
         ylabel('Y', 'FontSize', 12)
         title(sprintf('Problem #%d', i), 'FontSize', 12)
    end
    subplot(1,2,1)
    legend(handle', legNames)
end

%% Correlate parameters
% for i = 1:nProblems
%     ls(:,i) = clusterData{i}(:,1);
%     sv(:,i) = clusterData{i}(:,2);
%     nv(:,i) = clusterData{i}(:,3);
% end