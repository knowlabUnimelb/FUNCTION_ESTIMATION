function showSubjectData(fitsubject)
close all
datalocation = fullfile(pwd, 'Fits');
% datalocation = fullfile(pwd);
load(fullfile(datalocation, sprintf('gpAnalysis_full_s%d.mat', fitsubject)))

[m, midx] = max(fullfits.pM, [], 2);
table = [[1:8, 10:20, 22:26]', fullfits.pM, midx];

trialidx = [1:8 10:20 22:26]';

figure('WindowStyle', 'docked');
for i = 1:numel(trialidx)
    
    %% Set up data for fitting
    X     = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,4)), 4); % Observed X
    XSTAR = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,6)), 6); % Response X
    
    Y     = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,5)), 5); % Observed Y
    YSTAR = data(data(:,2) == trialidx(i,1) & ~isnan(data(:,7)), 7); % Response Y
    
    xyset       = (1:size([X, Y], 1))';
    [xy, xyord] = sortrows([X     Y    ], 1);
    xyset       = xyset(xyord);
    xystar      = sortrows([XSTAR YSTAR], 1);
    
    X = xy(:,1);         Y = xy(:,2);
    XSTAR = xystar(:,1); YSTAR = xystar(:,2);
    
    if ~isempty(X) && ~isempty(XSTAR) && numel(XSTAR) > 1
        respLocations = unique(floor(linspace(1, numel(XSTAR), 20)));
        xSTAR = XSTAR(respLocations);
        ySTAR = YSTAR(respLocations);
        
        subplot(5,6,i)
        plot(X, Y, ' ok', 'MarkerFaceColor', [0 0 0])
        hold on
%         plot(xSTAR, ySTAR, ' *k', 'MarkerSize', 5)
        plot(xSTAR, ySTAR, ' *r', 'LineWidth', 2)
        set(gca,'XLim', [-1 1], 'YLim', [-1 1], 'XTick', [], 'YTick', [])
        %         title(sprintf('%d - %d', trialidx(i), midx(i)))
        title(sprintf('%d', trialidx(i)))
    end
end