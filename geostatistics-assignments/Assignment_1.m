%% Assignment 1 GEO-SIGNAL ANALYSIS
%
% This is the script for the assignment 1 in Geo-signal analysis. It
% automatically loads the workspace variables I created, namely, station
% number, name and coordinates, as well as station numbers that have
% temperature measurements or not. Notice, that there are some stations
% that are missing from the latter case.
clear
clc

%% Definition of variable names:
%       LATnorth = y coordinate
%       LONeast = x coordinate
%       NAME = station name (city name)
%       STN = station number
%       STN_temp = station number that either has temperature measurement
%                  or has blank (no measurement)
%       temp = measured temperature
%
% At this point I have to erase all rows that correspond to NaNs in the temp
% column vector, as well as the rows that are missing according to
% STN_temp. I have to find the latter by comparing the station numbers
% (STN_temp and NAME column vectors)

load Ass_1_workspace

% Create a table such as: [station number, station name, y, x]
data = table(STN, NAME, LATnorth, LONeast);

% Create another table such as: [station number, temperature]
data_2 = table(STN_temp, temp);

% Sort rows of the tables in ascending order based on station number
data = sortrows(data, 1);
data_2 = sortrows(data_2, 1);

% I want to find which row of data corresponds to which row of data_2 and
% erase the remaining

% Find indices of elements that are not in
[values, idx_data, idx_data_2] = setxor(table2array(data(:, 1)), table2array(data_2(:, 1)));

% Erase the rows idx_data from data
data(idx_data, :) = [ ];

% Erase all the rows that correspond to NaNs in temp vector

% find the NaNs indices
idxNaN = find(isnan(table2array(data_2(:, 2))));

% find and store the station numbers that correspond to idxNaN so that I
% can erase them from data table as well
store = data_2(idxNaN, 1);

% Erase the rows idxNaN from data_2
data_2(idxNaN, :) = [ ];

% Erase the rows of stations that are stored to the store variable from
% data table. First I have to find in which positions those stations are
% place in the data table

[ ~, idxN, ~] = intersect(table2array(data(:,1)), table2array(store));

% Erase the rows
data(idxN, :) = [ ];

x = table2array(data(:,4)); % define x coordinates
y = table2array(data(:,3)); % define y coordinates
T = table2array(data_2(:,2)).*0.1; % define measured temperatures
%{
% Convert spherical to planar coordinates
x = deg2km(x)*10^3;
y = deg2km(y)*10^3;
%}
%% Export to xlsx the coordinates and the temperatures

% Multiply the temperatures by 10^-1 in order to take the true values
xlswrite('Temperature_coordinates.xlsx', [table2array(data(:,3:4)) 10^-1*table2array(data_2(:,2))])

%% Compute the Voronoi diagram

% Import the borders of the Netherlands and plot them
[lat lon] = borders('Netherlands');

voronoi(x, y)
title('Voronoi diagram of stations', 'Fontsize', 18)
xlabel('Longitude (decimal degrees)', 'Fontsize', 14);
ylabel('Latitude (decimal degrees)', 'Fontsize', 14);
axis square
hold on
plot(lon, lat, '-k', 'linewidth', 1.2);

%% Find average nearest neighbor distance
Near_N = zeros(length(x),1);
for i = 1:length(x)
    ap = distance(y(i), x(i), y, x); % Calculate distances among all points
    % Store nearest neighbor distance (the ap ~= 0 exludes the zero values
    % which are the distances between coincide points)
    Near_N(i) = min(ap(ap ~= 0));
end

% Find average nearest neighbor distance and convert it to kilometers
format bank
Near_N = deg2km(Near_N)*10^3 % convert to meters
NN_mean = mean(Near_N); % find the mean distance

fprintf('The average nearest neighbor distance is %.2f meters\n', NN_mean)

xlswrite('Nearest Neighbor distances.xlsx', Near_N)

%% Implement IDW correct


last = 500; % number of size of the grid I create below. This number will
% define the smoothness of the contours
multX = linspace(3.5, 7.2, last); % divide the x coordinates from minimum
X = ones(last, last); % define X as a matrix of ones
X = bsxfun(@times,multX, X); % multiply element by element the multX with the rows of X

multY = linspace(50.5, 53.5, last);
Y = ones(last, last); % define Y as a matrix of ones
Y = bsxfun(@times, multY, Y)'; % multiply element by element the multY with the rows of Y and find the transpose of it

temps = zeros(last, last); % preallocate the interpolated temperatures matrix
T = T'; % transform T to a row vector
p = 2; % power
% The following 2 loops calculate the temperature for each of the point of
% the grid whose coordinates are X and Y, based on the IDW interpolator
tic;
for i = 1:last
    for j = 1:last
        wd = 1./distance(Y(i, j), X(i, j), y, x).^p; % calculate the weight 1/d^p
        nom = T*wd; % linear combination of temperatures and wd
        denom = sum(wd); % sum of wd
        temps(i,j) = nom/denom; % calculate temperature
    end
end
timeIDW = toc;
% Plot the filled contours of the interpolated temperatures (temps matrix)
figure;
contourfm(Y, X, temps, 14)
colorbar
hold on

% Plot the station points
plot(x, y, '.r', 'markersize', 15)

% Import the borders of the Netherlands and plot them
plot(lon, lat, '-k', 'linewidth', 1.2);
xlim([3.5 7.25]); ylim([50.5 53.5]); axis square;
title('IDW interpolator of mean temperatures (C^o)', 'FontSize', 15);
xlabel('Longitude', 'fontsize', 14);
ylabel('Latitude', 'fontsize', 14);



%% Find best fit plane
% The plane is represented by T(x,y) = a*x + b*y + c + e <=>
% <=> T(x, y) = [x y 1][a; b; c] + e = <=> T(x, y) = A*[a; b; c] + e. By minimizing the error I get the
% least square estimates of the coefficients of the plain (a, b and c)

% Define design matrix A
A = [x y ones(34, 1)]; % A is 34x3

T = reshape(T,length(T),1);
% least square solution
estim = lscov(A, T);

% Compute the temperatures out of the estimated a, b and c coefficients
T_est = estim(1)*x + estim(2)*y + estim(3);

% Caclulate the detrended temperatures (residuals)
det_T = T - T_est;

% Visualize the original measured temperatures, the detrended temperatures
% and the best fit plain
last = 100;
% Create the grid of x and y
multX = linspace(3.5, 7.2, last); % divide the x coordinates from minimum
multY = linspace(50.5, 53.5, last);
[X, Y] = meshgrid(multX, multY);

% Compute the temperature for each point of the grid based on the plain
% function
bfitp = estim(1).*X + estim(2).*Y + estim(3);

figure;
% plot the measured temperatures
plot3(x, y, T, '.b', 'Markersize', 16);
hold on
% Plot the detrended temperatures (residuals)
plot3(x, y, det_T, '.g', 'Markersize', 16);

% plot the best-fit plain
surf(X, Y, bfitp);
shading interp
view(-140, 8)

xlabel('x coordinates', 'fontsize', 14);
ylabel('y coordinates', 'fontsize', 14);
zlabel('Temperature (C^o)', 'fontsize', 14);
leg = legend('Measured Temperatures', 'Detrended Temperatures (residuals)');
set(leg, 'fontsize', 15)
axis square
grid on


%% Experimental Variogram of original temperatures and detrended temperatures
[plot1, meanDist, meanSemiv] = ExperimentalVariogramIso(x, y, det_T, '.k', 'De', NN_mean)
[plot2, meanDist, meanSemiv] = ExperimentalVariogramIso(x, y, T, '.b', 'Me', NN_mean)

%% Biharmonic spline interpolation

last = 500;
multX = linspace(3.5, 7.2, last); % divide the x coordinates from minimum
multY = linspace(50.5, 53.5, last);

[X, Y] = meshgrid(multX, multY); % create the grid of X and Y coordinates

% Compute the interpolation values of the points of the grid. By 'v4' is
% declared that the biharmonic spline will be used
Z = griddata(x, y, T, X, Y, 'v4');

% Plot the interpolated surface
contourf(X, Y, Z, 10)
xlabel('x coordinates', 'Fontsize', 14)
ylabel('y coordinates', 'Fontsize', 14)
title('Biharmonic spline interpolation', 'Fontsize', 15)
