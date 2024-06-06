clear all% clc % close all%% Load hand rating data file data = dlmread('handParms.dat');% Massage raw data formatdata1 = data(data(:,1) ~= 2,:); % Remove the quadratic functionsdata1(data1(:,1)==4,1) = 2;     % Recode the simliarity functions at category 2data1(:,2:4) = data1(:,2:4);    % keep all parametersdata1(:,5) = ones(size(data1,1),1); % Add a column of ones%% Specify starting parameters and fit the linear planeraw_params = [0.35, -0.89, 0.37, 0.27, 0.79]; % [noise a1 a2 a3 b],start_params = norm_old_3dparams(raw_params); % Normalize starting parameters[final_params neglikelihood] = fit_3dGLC(start_params, data1, 7); % Fit linear discriminant planesubjlinbnd = [final_params(2:end)]; % Final parameters [a1 a2 a3 b],xyzaxes = [-4 4 -4 4 -7 0];range = [xyzaxes(1) xyzaxes(2); xyzaxes(3) xyzaxes(4); xyzaxes(5) xyzaxes(6)];%% Plot stimuli, optimal bound and other interesting boundsplot3dstim(data1, xyzaxes, 1); % Plot the parametersplot3dlinbnd(subjlinbnd,range(1:2,1:2),'r'); % Plot the boundaryxlabel('log Length Scale'); ylabel('log Signal Variance'); zlabel('log Noise Variance');grid on; hold onplot3(data(data(:,1) == 2, 2), data(data(:,1) == 2, 3), data(data(:,1) == 2, 4), 'or')predictions = double(lindiscrim3dvals(data1(:,2:4), subjlinbnd) < 0);predictions(lindiscrim3dvals(data1(:,2:4), subjlinbnd) >= 0) = 2;accuracy = percorr([data1(:,1:4), predictions], 3)