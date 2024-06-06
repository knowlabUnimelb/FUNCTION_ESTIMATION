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

    % Get the colors of all of the lines
    for i = 1:numel(h)
        c(i, :) = get(h(i), 'Color');
    end

    % Find the unique colors
    [u, idx] = unique(c, 'rows');

    % Delete black (which is the stimuli)
    idx(all(u == 0, 2), :) = [];
    c(all(c == 0, 2), :) = [];
    u(all(u == 0, 2), :) = [];

    % Count how many of each color
    for i = 1:size(u, 1)
        n(i) = sum(all(c == u(i,:), 2));
    end

    % sort them
    [n, nidx] = sort(n, 'descend');

    for i = 1:length(n)
        leglabels{i} = sprintf('Cluster %d (N = %d)', i, n(i));
    end
    legend(h(idx), leglabels, 'Location', 'best')

% end

set(gcf,'Position', [1.8        597.8        505.6        472.8])
set(gca, 'XTick',[], 'YTick', [])
xlabel('')
ylabel('') 
title('')

newfile = sprintf('newFig%02d.fig', dofig);
if exist(newfile, 'file') ~= 2
    savefig(newfile)
end