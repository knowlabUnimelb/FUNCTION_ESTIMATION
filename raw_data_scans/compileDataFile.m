clear all
% keep data
clc
close all

datafolder = 'C:\Documents and Settings\littled\My Documents\MATLAB\Analysis\Analysis SOCIAL FUNCTIONS\Gaussian Process Analysis\Data';
% datafolder = 'C:\Documents and Settings\littled\My Documents\My Dropbox\Work\Social Functions [gp analysis]\Data2';

% Need to scale data to axes [-1 1, -1 1]
data = dlmread('SocialFunctionData2.dat');

% Each of the 12 functions were presented once zoomed in and once zoomed
% out. Function 13 is a linearly increasing (or decreasing sine wave) and
% should be analysed separately
functionOrder = [1 11 3 4 6 3 5 11 13 8 7 12 2 9 8 10 12 10 4 2 13 7 1 6 9 5];

% zoom: 1 = in, 2 = out (F13: 1 = bottom left, 2 = bottom right)
      % 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
zoom = [1 1 1 1 1 2 1 2 2 2  2  1  2  2  1  1  2  2  2  1  1  1  2  2  1  2];

% in    out
xlims(:,:,1)  = [-1.2 1.2; -1.2 1.2; -1.2 1.2; -1.5 1.5; -1.2 1.2; -5.0  5.0; -1.2 1.2; -5.0 5.0; -3.0 3.0; -2.5 2.5; -2.5 2.5; -2.5 2.5; -6.5 16.];
ylims(:,:,1)  = [-1.6 1.6; -1.6 1.6; -1.0 1.0; -1.5 1.5; -0.8 1.6; -3.0 15.0; -2.2 1.2; -26. 6.0; -10. 20.; -16. 21.; -16. 21.; -7.0 15.;  1.5 30.];

xlims(:,:,2) = [-5.2 5.2; -3.2 3.2; -3.2 3.2; -3.5 3.5; -3.2 3.2; -8.0 8.0; -3.2 3.2; -10. 10.; -6.0 6.0; -7.5 7.5; -5.5 5.5; -7.5 7.5; -6.5 18.];
ylims(:,:,2) = [-5.6 5.6; -3.6 3.6; -3.0 3.0; -3.5 3.5; -3.5 4.1; -14. 26.0; -5.2 4.2; -45. 25.; -25. 35.; -30. 35.; -39. 39.; -25. 32.;  1.5 30.];

% Response minimum limits (max is always 1)
% rmins = dlmread('subjectMinimums.dat');
rmins = dlmread('subjectMinimums2.dat');

% data = dlmread('SocialFunctionData.dat');
subs = unique(data(:,1));
data(data(:,2) >= 26, 2) = data(data(:,2) >= 26, 2) - 1;
questions = unique(data(:,2));

funData = dlmread('functionList.dat');
%        FunIdx     FunClass        NData        Noise            X            Y

%% Sort the data file
% Columns:
% subject, trialcnt, responseCount, x, y, respX, respY, trialNumber, nDataPoints, GapLocation, nOutliers, Flipped?, Stimulus#
tic
for i = 1:numel(subs)
    newdata = [];
    for j = 1:numel(questions)
        d = data(data(:,1) == subs(i) & data(:,2) == questions(j), :);
        
        if ~isempty(d) %&& functionOrder(j) ~= 13
            fd = funData(funData(:,1) == functionOrder(j), :);            
            
            if functionOrder(j) == 13 && zoom(j) == 2
                fd(:,5) = fd(:,5) * -1 + 11.5;
            end
            
            % Rescale xy data
            xoffset = ((xlims(functionOrder(j),1,zoom(j)) - xlims(functionOrder(j),2,zoom(j)))/2 - xlims(functionOrder(j),1,zoom(j)));
            yoffset = ((ylims(functionOrder(j),1,zoom(j)) - ylims(functionOrder(j),2,zoom(j)))/2 - ylims(functionOrder(j),1,zoom(j)));
            xscale  = (xlims(functionOrder(j),2,zoom(j)) - xlims(functionOrder(j),1,zoom(j)))/2;
            yscale  = (ylims(functionOrder(j),2,zoom(j)) - ylims(functionOrder(j),1,zoom(j)))/2;
            
            x = (fd(:,5) + xoffset)./xscale;
            y = (fd(:,6) + yoffset)./yscale;
                       
            if rmins(i,1) == 0
                rxoffset = -.5; 
                ryoffset = -.5;
                rxscale = .5; 
                ryscale = .5;
            elseif rmins(i,1) == -1
                rxoffset = 0;
                ryoffset = 0; 
                rxscale = 1; 
                ryscale = 1; 
            end
            rx = (d(:,3) + rxoffset)./rxscale;
            ry = (d(:,4) + ryoffset)./ryscale;
            
            % Set up observed and response data
            if numel(x) <= size(d, 1)
                newdata = [newdata; ...
                    d(:,1), d(:,2),...
                    ones(size(d,1), 1),...
                    [[x y]; nan(size(d, 1) - size(fd, 1), 2)],...
                    [rx ry],...
                    d(:,1),...
                    repmat(fd(1,3), size(d, 1), 1),...
                    zeros(size(d,1),1),...
                    zeros(size(d,1),1),...
                    repmat(zoom(j), size(d,1), 1),...
                    repmat(functionOrder(j), size(d, 1), 1)];
            else
                newdata = [newdata; ...
                    d(:,1), d(:,2),...
                    ones(size(d,1), 1),...
                    [x y],...
                    [[rx ry]; nan(size(x,1) - size(d, 1), 2)],...
                    d(:,1),...
                    repmat(fd(1,3), size(d, 1), 1),...
                    zeros(size(d,1),1),...
                    zeros(size(d,1),1),...
                    repmat(zoom(j), size(d,1), 1),...
                    repmat(functionOrder(j), size(d, 1), 1)];
            end
        end
        disp(sprintf('Now processing question %d', j))
    end
    dlmwrite(fullfile(datafolder, sprintf('s%03d_socialFunctionData.dat', subs(i))), newdata, '\t')
    disp(sprintf('Now processing subject %d', i))
end
toc

%%
% ns = nSubplots(26);
% figure('WindowStyle', 'docked');
% for i = 1:26
%     subplot(ns(1), ns(2), i)
%     plot(newdata(newdata(:,2) == i,4), newdata(newdata(:,2) == i,5), ' ok'); hold on; 
%     plot(newdata(newdata(:,2) == i,6), newdata(newdata(:,2) == i,7), ' -k');
%     set(gca,'XLim', [-1 1], 'YLim', [-1 1])
% end