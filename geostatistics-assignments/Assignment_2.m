%% Assignment 2 GEO-SIGNAL ANALYSIS
%
% This is the script for the assignment 2 in Geo-signal analysis. It
% automatically loads the workspace variables I created, namely, station
% number, name and coordinates, as well as station numbers that have
% temperature measurements or not. Notice, that there are some stations
% that are missing from the latter case.
clear
clc
close all
% Definition of variable names:
%       x = x coordinates in degrees
%       y = y coordinates in degrees
%       T = Original temperatures
%       det_T = Detrended temperatures
%       NN_mean = Mean of the nearest neighbors

load Ass_2_variables

%+++++++ LEGEND ++++++++%
% plot1: Experimental variogram
% meanDist: Averages for every bin
% meanSemiv: Semivariance of the above meanDist
% matrix with the distances between observation points
[plot1, meanDist, meanSemiv] = ExperimentalVariogramIso(x, y, det_T, '.b', 'De', NN_mean);

% make the distances and the semviariances column vectors
meanDist = meanDist';
meanSemiv = meanSemiv';

%


%--------------%--------------Non-Linear Regression------------%------------%

%%%%%%%%%%%%%%%%%%%%%% Spherical %%%%%%%%%%%%%%%%%%%%%%%%

spher =@(c, dist) c(1).*(3.*dist./(2.*c(2)) - dist.^3./(2.*c(2).^3)).*(dist <= c(2)).*(dist > 0) + c(1).*(dist > c(2)) + 0.*(dist == 0);

% Initial values vector named ini. The first value is the sill, second is
% the range and third is the nugget
init = [var(det_T) 1.6*10^5];

% non-linear parameter estimation
[spher_hat, R_spher, J_spher, CovB_spher, MSE_spher] = nlinfit(meanDist, meanSemiv, spher, init);

dist_est_sph = feval(spher, spher_hat, meanDist);


% Compute values for the Spherical model in order to plot it
xx = 0:0.25:max(meanDist);
yy_sph = feval(spher, spher_hat, xx);

% Compute each point's standard deviation
er_spher = sqrt(diag(J_spher*CovB_spher*J_spher'));




%------------- Cross-Validation--------------%
    
% Create cell type. Each cell corresponds to a column, where each column
% consists of length(meanDist)-1 number of rows. In each column I have
% kicked out one cell. For example, from the 1st column I have kicked out
% the 1st cell, from the 2nd column the 2nd cell, from the 3rd column the
% 3rd cell etc. As a result, I will use the cell type in the
% cross-validation process

for i = 1:length(meanDist)
    crv_dist{i} = meanDist;
    crv_dist{i}(i) = [ ];
    crv_semiv{i} = meanSemiv;
    crv_semiv{i}(i) = [ ];
end


figure('Name', 'Spherical')
for i = 1:length(meanDist)
    % compute the estimates (sill, range) for the i-th cell kicked out
    crv_spher_hat(i, :) = nlinfit(crv_dist{i}, crv_semiv{i}, spher, init);
    % plot the predicted semivariances and the distances using the current
    % estimates (sill, range)
    plot(meanDist, feval(spher, crv_spher_hat(i,:), meanDist), '.', 'Markersize', 16)
    hold on
    % compute the predicted semivariance for the missing i-th cell
    crv_semiv{i} = feval(spher, crv_spher_hat(i,:), meanDist(i));
    % compute the difference (error) between the measured and predicted
    % semivariance
    crv_error_semiv(i) = crv_semiv{i} - meanSemiv(i);
end

% Compute Root Mean Squared Error (RMSE) of Cross-Validation
RMSE_spher = sqrt(mean(crv_error_semiv.^2));

% plot the above computed errors as error bars
%errorbar(meanDist, meanSemiv, crv_error_semiv, 'o')
% plot  the spherical model
p2 = plot(xx, yy_sph, '-b', 'Linewidth', 1.3);
% plot measurements
plot(meanDist, meanSemiv, '+k', 'Markersize', 13, 'Linewidth', 2);
xlabel('Distances (m)', 'Fontsize', 14)
ylabel('Semivariance', 'Fontsize', 14)
title('Spherical variogram model', 'Fontsize', 18)
grid on

message1 = sprintf('Nugget = 0\nSill = %.4f (C^o)^2\nRange = %.4f m', spher_hat(1), spher_hat(2))
text(10*10^4, 0.15, message1, 'Fontsize', 15, 'Backgroundcolor', 'c')


hold off


%%

%%%%%%%%%%%%%%%%%%%%%% Gaussian %%%%%%%%%%%%%%%%%%%%%%%%

gauss =@(c, dist) c(1).*(1 - exp(-(3.*dist).^2./c(2).^2));

% Initial values vector named ini. The first value is the sill, second is
% the range and third is the nugget
init = [var(det_T) 2*10^5];

% non-linear parameter estimation
[gauss_hat, R_gauss, J_gauss, CovB_gauss, MSE_gauss] = nlinfit(meanDist, meanSemiv, gauss, init);

dist_est_gauss = feval(gauss, gauss_hat, meanDist);

% Compute values for the Gaussian model in order to plot it
xx = 0:0.25:max(meanDist);
yy_gauss = feval(gauss, gauss_hat, xx);

% Compute each point's error
er_gauss = sqrt(diag(J_gauss*CovB_gauss*J_gauss'));



%------------- Cross-Validation--------------%

for i = 1:length(meanDist)
    crv_dist{i} = meanDist;
    crv_dist{i}(i) = [ ];
    crv_semiv{i} = meanSemiv;
    crv_semiv{i}(i) = [ ];
end


figure('Name', 'Gaussian')
for i = 1:length(meanDist)
    % compute the estimates (sill, range) for the i-th cell kicked out
    crv_gauss_hat(i,:) = nlinfit(crv_dist{i}, crv_semiv{i}, gauss, init);
    % plot the predicted semivariances and the distances using the current
    % estimates (sill, range)
    plot(meanDist, feval(gauss, crv_gauss_hat(i,:), meanDist), '.', 'Markersize', 16)
    hold on
    % compute the predicted semivariance for the missing i-th cell
    crv_semiv{i} = feval(gauss, crv_gauss_hat(i,:), meanDist(i));
    % compute the difference (error) between the measured and predicted
    % semivariance
    crv_error_semiv(i) = crv_semiv{i} - meanSemiv(i);
end

% Compute Root Mean Squared Error (RMSE) of Cross-Validation
RMSE_gauss = sqrt(mean(crv_error_semiv.^2));

% plot the above computed errors as error bars
%errorbar(meanDist, meanSemiv, crv_error_semiv, 'o')
% plot  the gaussian model
p2 = plot(xx, yy_gauss, '-g', 'Linewidth', 1.3);
% plot measurements
plot(meanDist, meanSemiv, '+k', 'Markersize', 13, 'Linewidth', 2);
xlabel('Distances (m)', 'Fontsize', 14)
ylabel('Semivariance', 'Fontsize', 14)
title('Gaussian variogram model', 'Fontsize', 18)
grid on

message2 = sprintf('Nugget = 0\nSill = %.4f (C^o)^2\nRange = %.4f m', gauss_hat(1), gauss_hat(2))
text(10*10^4, 0.15, message2, 'Fontsize', 15, 'Backgroundcolor', 'c')


hold off

%

%%%%%%%%%%%%%%%%%%%%%% Exponential %%%%%%%%%%%%%%%%%%%%%%%%

expon =@(c, dist) c(1).*(1 - exp(-3.*dist./c(2)));

% Initial values vector named ini. The first value is the sill, second is
% the range and third is the nugget
init = [var(det_T) 1.5*10^5];

% non-linear parameter estimation
[expon_hat, R_expon, J_expon, CovB_expon, MSE_expon] = nlinfit(meanDist, meanSemiv, expon, init);

dist_est_expon = feval(expon, expon_hat, meanDist);

% Compute values for the Exponential model in order to plot it
xx = 0:0.25:max(meanDist);
yy_expon = feval(expon, expon_hat, xx);

% Compute each point's error
er_expon = sqrt(diag(J_expon*CovB_expon*J_expon'));



%------------- Cross-Validation--------------%

for i = 1:length(meanDist)
    crv_dist{i} = meanDist;
    crv_dist{i}(i) = [ ];
    crv_semiv{i} = meanSemiv;
    crv_semiv{i}(i) = [ ];
end


figure('Name', 'Exponential')
for i = 1:length(meanDist)
    % compute the estimates (sill, range) for the i-th cell kicked out
    crv_expon_hat(i,:) = nlinfit(crv_dist{i}, crv_semiv{i}, expon, init);
    % plot the predicted semivariances and the distances using the current
    % estimates (sill, range)
    plot(meanDist, feval(expon, crv_expon_hat(i,:), meanDist), '.', 'Markersize', 16)
    hold on
    % compute the predicted semivariance for the missing i-th cell
    crv_semiv{i} = feval(expon, crv_expon_hat(i,:), meanDist(i));
    % compute the difference (error) between the measured and predicted
    % semivariance
    crv_error_semiv(i) = crv_semiv{i} - meanSemiv(i);
end

% Compute Root Mean Squared Error (RMSE) of Cross-Validation
RMSE_expon = sqrt(mean(crv_error_semiv.^2));

% plot the above computed errors as error bars
%errorbar(meanDist, meanSemiv, crv_error_semiv, 'o')
% plot  the exponential model
p2 = plot(xx, yy_expon, '-m', 'Linewidth', 1.3);
% plot measurements
plot(meanDist, meanSemiv, '+k', 'Markersize', 13, 'Linewidth', 2);
xlabel('Distances (m)', 'Fontsize', 14)
ylabel('Semivariance', 'Fontsize', 14)
title('Exponential variogram model', 'Fontsize', 18)
grid on

message3 = sprintf('Nugget = 0\nSill = %.4f (C^o)^2\nRange = %.4f m', expon_hat(1), expon_hat(2))
text(10*10^4, 0.15, message3, 'Fontsize', 15, 'Backgroundcolor', 'c')


hold off

%
% plot models
figure('Name', 'All models');
plot(meanDist, meanSemiv, '+k', 'Markersize', 13, 'Linewidth', 3);
hold on
% spherical
p2 = plot(xx, yy_sph, '-b', 'Linewidth', 1.1);
% gauss
p3 = plot(xx, yy_gauss, '-g', 'Linewidth', 1.1);
% gauss
p4 = plot(xx, yy_expon, '-m', 'Linewidth', 1.1);
grid on

% plot error bars for every model
errorbar(meanDist, meanSemiv, er_spher, 'ob')
errorbar(meanDist, meanSemiv, er_gauss, 'og')
errorbar(meanDist, meanSemiv, er_expon, 'om')

xlabel('Distances', 'Fontsize', 14)
ylabel('Dissimilarity', 'Fontsize', 14)
grid on
legend([p2 p3 p4], {'Spherical', 'Gaussian', 'Exponential'})

% Plot normal QQ plots of the residuals
figure('Name', 'Residuals');
subplot(2,3,1)
qqplot(R_spher), title('Spherical Residuals', 'Fontsize', 14)
grid on, axis square
subplot(2,3,4)
histfit(R_spher, 4), grid on, axis square
subplot(2,3,2)
qqplot(R_gauss), title('Gaussian Residuals', 'Fontsize', 14)
grid on, axis square
subplot(2,3,5)
histfit(R_gauss, 4), grid on, axis square
subplot(2,3,3)
qqplot(R_expon), title('Exponential Residuals', 'Fontsize', 14)
grid on, axis square
subplot(2,3,6)
histfit(R_expon, 4), grid on, axis square


% Present results

MODEL = [{'Spherical', 'Gaussian', 'Exponential'}]';
MSE_nlinfit = [MSE_spher MSE_gauss MSE_expon]';
Norm_Residuals_nlinfit = [norm(R_spher) norm(R_gauss) norm(R_expon)]';
RMSE_CrossValidation = [RMSE_spher RMSE_gauss RMSE_expon]';
Results = table(MODEL, MSE_nlinfit, Norm_Residuals_nlinfit, RMSE_CrossValidation)

%------ Clear JUNK --------%

clear message1 message2 message3 p1 p2 p3 p4

%% Question 2
%close all
%clc
%------------------- Spherical Covariance Function ----------------%

% Define the spherical covariance function
cov_spher =@(c, dist) c(1).*(1 - 3.*dist./(2.*c(2)) + dist.^3./(2.*c(2).^3)).*(dist>=0).*(dist<=c(2)) + 0.*(dist>c(2));

% Evaluate the spherical covariance function
cov_yy_sph = feval(cov_spher, spher_hat, xx);

% Plot
figure('Name', 'Spherical Covariance Function');
plot(xx, cov_yy_sph, '-b', 'Linewidth', 1.3);
xlabel('Distance (m)', 'Fontsize', 14)
ylabel('Covariance (C^o)^2', 'Fontsize', 14)
title('Spherical Covariance Function', 'Fontsize', 18);
grid on

set(gca, 'ytick', [0.05:0.05:0.25 spher_hat(1)]);
set(gca, 'xtick', [0:2e+4:4e+4 spher_hat(2) 6e+4:2e+4:18e+4], 'xticklabelrotation', 45);

text(0.7e+4, 0.275, '\leftarrowSill', 'Fontsize', 18, 'Backgroundcolor', 'c')
text(5.6e+4, 0.02, '\downarrowRange', 'Fontsize', 18, 'Backgroundcolor', 'c')

%% Question 3 - Find the OK weights
close all

%------ Compute the Ordinary Kriging weights and estimations -------%

[OKw, T_est] = ordinarykriging(mean(x), mean(y), x, y, T, spher, spher_hat, 1);

% Erase the Lagrange multiplier
OKw(end) = [ ];

fprintf('The sum of the OK weights is %d\n', sum(OKw))

% Plot the mean(x) mean(y) point with its weights
figure('Name', 'Weights of the mean coordinates');
scatter(x, y, 100, OKw, 'filled');
colorbar
hold on
p1 = plot(mean(x), mean(y), 'ro', 'Linewidth', 2.2, 'Markersize', 10);
xlabel('Longitude', 'Fontsize', 14)
ylabel('Latitude', 'Fontsize', 14)
title('Weights w.r.t. mean x and y', 'Fontsize', 18)
grid on
axis square

% Plot borders of Netherlands
% Import the borders of the Netherlands and plot them
[lat lon] = borders('Netherlands');

plot(lon, lat, '-k', 'linewidth', 1.2);

%----- Find the three closes stations -----%

d = distance(y, x, mean(y), mean(x));
% Assign number to each station
d = [[1:length(x)]' d];

% Sort in ascending order
d = sortrows(d, 2);

% Plot the three closes points
p2 = plot(x(d(1:3,1)), y(d(1:3,1)), 'squarer', 'Linewidth', 2.2);
leg = legend([p1 p2], {'Mean Center', '3 Closest Points'}, 'Location', 'Northwest');
set(leg, 'Fontsize', 14)

% Plot the Temperatures of the 3 closest points
message1 = sprintf(['T_1 = ' num2str(T(d(1,1))) '^oC\nT_2 = ' num2str(T(d(2,1))) '^oC\nT_3 = ' num2str(T(d(3,1))) '^oC\nT_m = ' num2str(T_est) '^oC']);
text(3.25, 53, message1, 'Backgroundcolor', 'c', 'Fontsize', 14)

%% Question 4 - 
close all

%------ Compute the Ordinary Kriging weights and estimations -------%

[OKw, T_est, variK] = ordinarykriging(x, y, x, y, T, spher, spher_hat, 1);

% Plot estimated and measured temperatures
figure('Name', 'Temperatures');
plot(1:length(T), T, '-b', 'Linewidth', 1.2);
hold on
plot(1:length(T_est), T_est, '--r', 'Linewidth', 1.2)
xlabel('Number of point', 'Fontsize', 14)
ylabel('Temperature (^oC)', 'Fontsize', 14)
grid on
leg = legend('Measured', 'Estimated');
set(leg, 'Fontsize', 14)

% Compute error
er = mean(T-T_est)


%% Question 5 - Estimate Temperature and Variance at a grid (of IDW)

last = 500; % number of size of the grid I create below. This number will
            % define the smoothness of the contours
multX = linspace(3.5, 7.2, last); % divide the x coordinates from minimum
multY = linspace(50.5, 53.5, last);

[X, Y] = meshgrid(multX, multY);

%------ Compute the Ordinary Kriging weights and estimations -------%

[OKw, T_est, variK, timeOK] = ordinarykriging(X, Y, x, y, T, spher, spher_hat, 1);

% Plot the interpolated map
figure('Name', 'Ordinary Kriging Temperatures');
contourf(X, Y, T_est,14);
colorbar
hold on
plot(lon, lat, '-k', 'linewidth', 1.2);
plot(x, y, '.r', 'Markersize' ,15);
xlabel('Longitude', 'Fontsize', 14)
ylabel('Latitude', 'Fontsize', 14)
title('Ordinary Kriging Temperatures', 'Fontsize', 18)
axis square
xlim([3.5 7.2])

% Plot the variances map
figure('Name', 'Ordinary Kriging Variances')
subplot(1,2,1)
contourf(X, Y, variK, 15);
hold on;
plot(x, y, '.r', 'Markersize', 15);
colorbar
axis square

subplot(1,2,2)
image(variK, 'CDataMapping', 'scaled')
axis square
colorbar

%% Question 6 - Differences between same grid (500x500) of IDW and OKriging

% Load estimated temperatures with IDW and time it take to execute
load Time_TempsIDW.mat

IDWOKdiff = abs(temps - T_est);

figure('Name', 'Variance - Difference IDWOK');
plot(variK, IDWOKdiff, '.r');
xlabel('Variances', 'Fontsize', 14);
ylabel('Differences', 'Fontsize', 14);

figure('Name', 'Differences IDW - OK');

subplot(1,2,1)
contourf(X, Y, IDWOKdiff, 10);
colorbar
grid on
axis square
title('Differences IDW - OK', 'Fontsize', 18)
xlabel('Longitude', 'Fontsize', 14);
ylabel('Latitude', 'Fontsize', 14);
hold on
plot(x, y, '.r', 'Markersize', 13)

subplot(1,2,2)
image(IDWOKdiff, 'CDataMapping', 'scaled')
colorbar
title('Differences IDW - OK', 'Fontsize', 18);
axis square

%% Question 7 - Indicator

% Define indicator heights between min and max of T
ind_heights = linspace(min(T), max(T), 10);

% Calculate percentages
perc_indic = sum(indic)./length(T)*100;

% Plot the CDF of temperatures
figure('Name', 'Empirical Cumulative Distribution Function');
cdfplot(T)

%----------- Find Threshold Temperatures for 60% and 90% ------------%

% Extract the Temperatures and their corresponding probabilities
[prob T_crits] = ecdf(T)

% Find Temperature for 60% probability using linear interpolation
T_crits60 = interp1(prob, T_crits, 0.6)

% Find Tepmerature for 90% probability using linear interpolation
T_crits90 = interp1(prob, T_crits, 0.9)

%%
% Convert the observation to INDICATOR values

%~~~~~~~ 60% ~~~~~~~~~%

% Pre-allocation of indic
indic_60 = zeros(size(T));

% Determine Indicator height for 60%
indic_60 = 1.*(T <= T_crits60) + 0.*(T > T_crits60);

%~~~~~~~~ 90% ~~~~~~~~~%

% Pre-allocation of indic
indic_90 = zeros(size(T));

% Determine Indicator height for 90%
indic_90 = 1.*(T <= T_crits90) + 0.*(T > T_crits90);


%% Compute Indicator experimental variogram

[X1, X2] = meshgrid(x);
[Y1, Y2] = meshgrid(y);
[T1, T2] = meshgrid(indic_60);

D = deg2km(distance(Y2, X2, Y1, X1))*10^3; % Calculate distances among stations and convert them to meters
G = 0.5*(T1 - T2).^2; % Calculate dissimilarity (semivariance)

index = 1:length(indic_60); % vector with the number of points
% create two matrices out of index vector where the one is rotated 90 degrees
[L1 L2] = meshgrid(index);
% Build the logical matrix
I = L1 > L2;

% Plot the variogram
figure;
plot(D(I), G(I), '.r', 'Markersize', 12); % plot the true values
xlabel('Distance (m)', 'Fontsize', 14);
ylabel('Dissimilarity', 'Fontsize', 14);

%%

heights = linspace(min(T), max(T), 15);

for i = 1:length(heights)
    in(:,i) = 1.*(T<=heights(i)) + 0.*(T>heights(i));
end

s = sum(in)./34;