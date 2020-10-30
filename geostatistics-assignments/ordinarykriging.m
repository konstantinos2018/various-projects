function [ OKw, T_est, variK, timeOK ] = ordinarykriging( x_est, y_est, x, y, value, model, params )
% Function ordinarykriging computes the OK weights that are going to be
% used in the Ordinary Kriging interpolation technique.
%
% ___INPUTS___
%
% semiv: Semivariance as computed by a variogram model (either matrix or array)
% x_est: Array of x coordinates of the points where the estimation
%        will take place
% y_est: Array or matrix of y coordinates of the points where the estimation
%        will take place
% x: Array of x coordinates of the observation points
% y: Array of y coordinates of the observation points
% model: function_handle with the model function (e.g. @(x) 1./x.^2)
% params: vector consists of the model parameters with the order [sill
%         range]
% grd_rand_points = either 'GridPoints' (0) or 'RandomPoints' (1). Estimate values
% either at gridded points or random points
%
%____OUTPUTS____
%
% OKw: Ordinary Kriging weights
% T_est: Estimated values


%===================== BEGIN ======================%

% find the distances among the observation points
[x1, x2] = meshgrid(x);
[y1, y2] = meshgrid(y);

dist_diff = distance(y2, x2, y1, x1);

% Plot the distance matrix
figure('Name', 'Distances of Observation points')
image(dist_diff, 'cdatamapping', 'scaled')
colorbar
title('Distances of observation points')

%------------- DESIGN MATRIX -----------%

% Compute semivariances from variogram function of OBSERVATION points (This
% is the DESIGN matrix)
sem = feval(model, params, dist_diff);

% Plot the above semivariances
figure('Name', 'Semivariances of Observation points')
image(sem, 'cdatamapping', 'scaled')
colorbar
title('Semivariances of observation points')

% Finalize the DESIGN matrix by concatenating the ones to the bottom and
% the far-right and a zero at the sem(end,end)

sem = [sem; ones(1, length(x))];
sem = [sem ones(length(x)+1, 1)];
sem(length(x)+1, length(x)+1) = 0;

figure('Name', 'Semivariances with Lagrange');
image(sem, 'CDataMapping', 'scaled')
colorbar
title('Semivariances FULL')

% Pre-allocate T_est and variK
T_est = zeros(size(x_est));
variK = zeros(size(x_est));

if isvector(x_est) & isvector(y_est)

    %--------------- INTERPOLATION ----------------%

    for i = 1:length(x_est)
        % compute distances between the observation points and the estimation
        % point
        dist_RH = distance(y_est(i), x_est(i), y, x);

        % Compute the semivariance for the above distances from variogram
        % function
        RH = feval(model, params, dist_RH);

        % Concatenate the number one that is assigned to the last cell of
        % the array, which is the sum of the weights (sum of weights equal
        % to one)
        RH(length(x)+1) = 1;

        % Compute the weights
        OKw = sem\RH;

        % Compute estimations which is the linear combination of
        % detrended temperatures and the weights
        T_est(i) = OKw(1:length(x))'*value;

        % Compute variance
        variK(i) = OKw(1:length(x))'*RH(1:length(x)) + OKw(length(x)+1);
    
    end
    
elseif ~isvector(x_est) & ~isvector(y_est)
    
    %--------------- INTERPOLATION ----------------%
    
    % store size of grid x_est
    s = size(x_est);
    
    % Reshape x_est and y_est and make them vectors
    X = reshape(x_est, [], 1);
    Y = reshape(y_est, [], 1);
    tic;
    for i = 1:length(X)
        % compute distances between the observation points and the estimation
        % point
        dist_RH = distance(Y(i), X(i), y, x);

        % Compute the semivariance for the above distances from variogram
        % function
        RH = feval(model, params, dist_RH);

        % Concatenate the number one that is assigned to the last cell of
        % the array, which is the sum of the weights (sum of weights equal
        % to one)
        RH(length(x)+1) = 1;

        % Compute the weights
        OKw = sem\RH;

        % Compute estimations which is the linear combination of
        % detrended temperatures and the weights
        T_est(i) = OKw(1:length(x))'*value;

        % Compute variance
        variK(i) = OKw(1:length(x))'*RH(1:length(x)) + OKw(length(x)+1);
    
    end
    timeOK = toc;
    % reshape the T_est and varik in the size of x_est
    T_est = reshape(T_est, s);
    variK = reshape(variK, s);
    
end     