function [ var_plot, DE, GE ] = ExperimentalVariogramIso(x_coord, y_coord, value, color, type, near_nei)
%%% Experimental Variogram function
%   Legend of variables:
%       x_coord = vector of the x coordinates
%       y_coord = vector of the y coordinates
%       value = vector of the measured quantity (temperature)
%       color = color of the points of the variogram
%       type = type of data ('Me' or 'De' for measured and detrended)
    
    % X2 (Y2 and T2) is X1 (Y1 and T1) rotated by 90 degrees
    [X1, X2] = meshgrid(x_coord);
    [Y1, Y2] = meshgrid(y_coord);
    [T1, T2] = meshgrid(value);
    
    D = deg2km(distance(Y2, X2, Y1, X1))*10^3; % Calculate distances among stations and convert them to meters
    G = 0.5*(T1 - T2).^2; % Calculate dissimilarity (semivariance)
    
    % Process that will take into account only the lower diagonal (or it
    % could take the upper diagonal) matrices of D and G by building a
    % logical matrix in the following way:
    index = 1:length(value); % vector with the number of points
    % create two matrices out of index vector where the one is rotated 90 degrees
    [L1 L2] = meshgrid(index);
    % Build the logical matrix
    I = L1 > L2;
    
    % Plot the variogram
    figure;
    plot(D(I), G(I), color, 'Markersize', 12); % plot the true values
    xlabel('Distance (m)', 'Fontsize', 14);
    ylabel('Dissimilarity (C^o)^2', 'Fontsize', 14);
    
    if type == 'Me'
        title('Measured mean temperature variogram', 'Fontsize', 15);
    elseif type == 'De'
        title('Detrended mean temperature variogram', 'Fontsize', 15);
    else
        error('There is no such an option for the type of data')
    end
    grid on
    hold on
    
    % Group the distances into bins and take the average of each group
    % 
    % Set the lag equal to the mean nearest neighbor distance
    lag = near_nei;
    lag = 18000;
    
    % As the maximum distance of the variogram, as a rule of thumb, I take
    % the half of the maximum of the distances
    hmd = max(D(:))/1.6;    
    % Now, the number of the groups (lags) is the maximum distance divided
    % by the lag magnitude 
    n_lags = floor(hmd/lag);
    fprintf('Number of lags = %d\n', n_lags)
    LAGS = ceil(D/lag);
    
    % compute the mean value in each group
    for i = 1:n_lags
        in_group = (LAGS == i); % logical index that shows to which group each pair of distance belongs to
        % Finds the mean of the distances that belong to lag i
        DE(i) = mean(mean(D(in_group)));
        % finds the number of pairs that belong to the group i
        PN(i) = sum(sum(in_group == 1))/2;
        GE(i) = mean(mean(G(in_group))); % compute mean semivariance of group i
    end
    
    % Plot the mean of distances of each group against the corresponding
    % mean semivariance
    var_plot = plot(DE, GE, '+r', 'Markersize', 14, 'linewidth', 3);
    xlim([0 hmd]);
    % Plots the horizontal line that crossed the sill value and the
    % variance
    vT = var(value);
    xx = [0 max(DE)];
    yy = [vT vT];
    plot(xx, yy, '--g', 'linewidth', 2)
end
