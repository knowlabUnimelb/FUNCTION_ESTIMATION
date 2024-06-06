function showAllSubjectsData(fititem)
% fititem = 10;

close all
datalocation = fullfile(pwd, 'fits');
subjectFileNames = dir(fullfile(datalocation, 'gpAnalysis_full_s*.mat'));

rX = nan(40, numel(subjectFileNames));
rY = nan(40, numel(subjectFileNames));

tic
for subjectIdx = 1:numel(subjectFileNames)
    load(fullfile(datalocation, subjectFileNames(subjectIdx).name))
    
    if subjectIdx == 1;
        oX = fullfits.X{fititem};
        oY = fullfits.Y{fititem};
    end
    
    if numel(fullfits.xSTAR{fititem}) >= 40
        rX(:,subjectIdx) =  fullfits.xSTAR{fititem}(1:40);
    else
        rX(1:numel(fullfits.xSTAR{fititem}), subjectIdx) = fullfits.xSTAR{fititem};
    end
    if numel(fullfits.ySTAR{fititem}) >= 40
        rY(:,subjectIdx) =  fullfits.ySTAR{fititem}(1:40);
    else
        rY(1:numel(fullfits.ySTAR{fititem}), subjectIdx) = fullfits.ySTAR{fititem};
    end
    [m, modelidx(subjectIdx)] = max(fullfits.pM(fititem, :));
end
displayTime(toc);

colors = [1 0 0; .75 .25 0; .5 .5 1; .25 .75 1; 0 1 0]; 

figure('WindowStyle', 'docked');
modelCounts = aggregate([ones(numel(modelidx), 1), modelidx'], 2, 1, @count);
smc = sortrows(modelCounts, -2);
mnames = {'Linear', 'Quadratic', 'Cubic', 'Quartic', 'Similarity'};
for i = 1:size(smc, 1)
    h = plot(rX(:,modelidx == smc(i,1)), rY(:,modelidx == smc(i,1)), '-', 'Color', colors(smc(i,1),:));
    lh(i) = h(1);
    legnames{i} = sprintf('%s - N = %d', mnames{smc(i,1)}, smc(i,2));
    hold on
end
plot(oX, oY, ' ok', 'MarkerFaceColor', [0 0 0], 'MarkerSize', 10)
set(gca, 'XLim', [-1 1], 'YLim', [-1 1], 'XTick', [], 'YTick', [])
legend(lh, legnames)
