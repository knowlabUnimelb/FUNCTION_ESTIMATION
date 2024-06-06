clear all
clc

load parameterClusters
keep subIdx clusterData

load handRatings

annieIdx = annie(:,1); annie(:,1) = []; jerry(:,1) = []; 
annie(annie == 3) = 4; jerry(jerry == 3) = 4; 
annie(annie == 5 | annie == 6) = nan;
jerry(jerry == 5 | jerry == 6) = nan;

nanSubs = unique([find(all(isnan(annie), 2)); find(all(isnan(jerry), 2))]);
ratedSubs = annieIdx;
ratedSubs(ismember(annieIdx, nanSubs)) = [];
annie(ismember(annieIdx, nanSubs),:) = [];

for i = 1:size(annie,2)
    for j = 1:size(annie, 1)
        if ~isempty(find(subIdx == ratedSubs(j)))
        ls(j,i) = clusterData{i}(find(subIdx == ratedSubs(j)), 1);
        sv(j,i) = clusterData{i}(find(subIdx == ratedSubs(j)), 2);
        nv(j,i) = clusterData{i}(find(subIdx == ratedSubs(j)), 3);
        else
        ls(j,i) = nan;
        sv(j,i) = nan;
        nv(j,i) = nan;
        end
    end
end

%% Scatter all parameters from all functions
subNums = repmat((1:size(annie,1))', 1, size(annie,2));
funNums = repmat(1:size(annie,2), size(annie,1), 1);
bigMat = [subNums(:), funNums(:), annie(:), ls(:), sv(:), nv(:)];
bigMat(any(isnan(bigMat), 2), :)=[];

markers = {' xb', ' og', ' sr', ' ^c', ' +m', '*k'};
lineColors = {'b', 'g', 'r', 'c', 'm', 'k'};
clusterNames = unique(bigMat(:,3));
for i = 1:numel(clusterNames); 
    subplot(2,2,1)
    h1 = plot3(bigMat(bigMat(:,3) == clusterNames(i), 4),...
               bigMat(bigMat(:,3) == clusterNames(i), 5),...
               bigMat(bigMat(:,3) == clusterNames(i), 6),...
               markers{i});
    xlabel('log Length Scale', 'FontSize', 12)
    ylabel('log Signal Variance', 'FontSize', 12)
    zlabel('log Noise Variance', 'FontSize', 12)
    hold on
    box on

    subplot(2,2,2)
    h2 = plot(bigMat(bigMat(:,3) == clusterNames(i), 4), bigMat(bigMat(:,3) == clusterNames(i), 5), markers{i});
    hold on
    xlabel('log Length Scale', 'FontSize', 12)
    ylabel('log Signal Variance', 'FontSize', 12)
    
    subplot(2,2,3)
    h3 = plot(bigMat(bigMat(:,3) == clusterNames(i), 4), bigMat(bigMat(:,3) == clusterNames(i), 6), markers{i});
        hold on
    xlabel('log Length Scale', 'FontSize', 12)
    ylabel('log Noise Variance', 'FontSize', 12)
    
    subplot(2,2,4)
    h4 = plot(bigMat(bigMat(:,3) == clusterNames(i), 5), bigMat(bigMat(:,3) == clusterNames(i), 6), markers{i});
        hold on
    xlabel('log Signal Variance', 'FontSize', 12)
    ylabel('log Noise Variance', 'FontSize', 12)
end