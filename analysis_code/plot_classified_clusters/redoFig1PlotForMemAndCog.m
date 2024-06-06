% function redoPlotForMemAndCog(probnumber)
close all
clear all
clc

dofig = 24; 

probIndex = [21,  6,  1, 3, 24, 10, 7, 20, 13,  8, 23,  2,...
             12, 18, 19, 4, 22,  9, 5, 14, 17, 16, 15, 11];


probnumber = probIndex(dofig); 

% Open specified figure (wll contain two supblots (cluster and hand drawn
% functions)
uiopen(fullfile(pwd, sprintf('Prob%02d.fig',probnumber)), 1)

% Doc the window
% set(gcf, 'WindowStyle', 'docked')

% Find and delete left hand subplot (clusters)
hf = get(gcf,'Children');

switch probnumber
    case {6, 24, 10, 7, 13, 8, 23, 12, 18, 19, 4, 22, 5, 14, 11, 16}
        delete(hf(3))
    case {21, 3, 20}
        delete(hf(1))
        delete(hf(4))
        hf(1) = hf(2);
    case {1, 2, 9, 15, 17}
        delete(hf(1))
        delete(hf(4))
        hf(1) = hf(2);
    otherwise
        hf
        return
end

% Reposition functions
set(hf(1), 'Position', [0.17778      0.17597      0.72722      0.74903])

% change font size
set(hf(1), 'FontSize', 16)

% If there isn't a legend already, create one
% if isempty(get(gca, 'Legend'))
    h = get(hf(1), 'Children'); % Find all of the lines

    if dofig == 16
        error('this one broken. do by hand.')
        return
    else
    figure; 
    plot(get(h(1), 'XData'), get(h(1), 'YData'), 'Color', get(h(1), 'Color'),...
        'LineStyle', 'none', 'Marker', 'o', 'MarkerFaceColor', 'k', 'MarkerSize', 10); 
    set(gca,'XLim', [-1 1], 'YLim', [-1 1])
    end

% end

set(gcf,'Position', [1.8        597.8        505.6        472.8])
set(gca, 'XTick',[], 'YTick', [])
xlabel('')
ylabel('') 
title('')

newfile = sprintf('newStimulus_Fig%02d.fig', dofig);
if exist(newfile, 'file') ~= 2
    savefig(newfile)
end