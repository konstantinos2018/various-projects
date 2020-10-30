%%

clear all
close all
clc

%%

%%%% Load data %%%%
data = load('Torfajokull_defo_env.txt');
deform = data(:,5); % deformation data to deform array
procrd = data(:,3:4); % projected coordinates to 10000x2 matrix (in meters)
data(:, 3:5) = [ ];
%%

%%%%%%%%%% Questions 1 %%%%%%%%%%

% find indeces of samples that show uplift/subsidence
uplift = find(deform > 0); % find cells with positive deformation
subsid = find(deform < 0); % find cells with negative deformation

%%%% Spatial visualisation of uplift and subsidence %%%%
figure('Name', 'x-y-deformation scatter', 'numbertitle', 'off');

ax1 = subplot(1,2,1);
scatter(procrd(uplift,1), procrd(uplift,2), 15, deform(uplift), 'filled');
xlabel('X coordinates (m)', 'fontsize', 13);
ylabel('Y coordinates (m)', 'fontsize', 13);
title('Uplift spatial pattern', 'fontsize', 15);
grid on;
axis square;
c1 = colorbar; % display colorbar
colormap(ax1, parula); % define the colors of the bar
c1.Label.String = 'm/y'; % label the colorbar
c1.Label.FontSize = 13; % fontsize of label


ax2 = subplot(1,2,2);
scatter(procrd(subsid,1), procrd(subsid,2), 15, deform(subsid), 'filled');
xlabel('X coordinates (m)', 'fontsize', 13);
ylabel('Y coordinates (m)', 'fontsize', 13);
title('Subsidence spatial pattern', 'fontsize', 15);
grid on;
axis square;
c2 = colorbar;
colormap(ax2, flipud(parula)); % define inverted colors of the bar
c2.Label.String = 'm/y';
c2.Label.FontSize = 13;

%%

%%%%%%%%%% Question 3 %%%%%%%%%%
figure
scatter(procrd(subsid,1), procrd(subsid,2), 15, deform(subsid), 'filled')
xlabel('X coordinates (m)', 'fontsize', 13);
ylabel('Y coordinates (m)', 'fontsize', 13);
title('Study area of subsidence', 'fontsize', 16);
grid on
axis square
hold on

%%%% Initial value of DEPTH %%%%

%%%%%%%%%%%%%%% First assumption about depth %%%%%%%%%%%%%%%%%%

% WRONG initial value and circle of the WHOLE study area

minx = min(procrd(subsid, 1)); % far West boundary of study area
maxx = max(procrd(subsid, 1)); % far East boundary of study area
miny = min(procrd(subsid, 2)); % far South boundary of study area
maxy = max(procrd(subsid, 2)); % far North boundary of study area

radius = (maxx - minx)/2; % radius of circle
theta = [0:pi/40:2*pi]; % angle of the circle (0 to 2*pi)

% midx and midy indicate the displacement of the circle in the x, y
% coordinates from the origin
midx = (maxx + minx)/2; % displacement in x
midy = (maxy + miny)/2; % displacement in y

xcir = radius*cos(theta) + midx; % x coordinates of circle
ycir = radius*sin(theta) + midy; % y coordinates of circle

p1 = plot(xcir, ycir, 'k', 'linewidth', 1.8); % plot circle


%%%%%%%%%%%%%%%%%%%% Second assumption about depth %%%%%%%%%%%%%%%%%%%%

% Find and plot circle that almost surrounds the subsided study area

subsid2 = find(deform <= -8e-3); % threshold defined as -8e-3

xdisp = abs(deform(subsid2))'*procrd(subsid2,1)/sum(abs(deform(subsid2)));
ydisp = abs(deform(subsid2))'*procrd(subsid2,2)/sum(abs(deform(subsid2)));

radius = 5000; % radius of circle
theta = [0:pi/40:2*pi]; % angle of the circle (0 to 2*pi)

% midx and midy indicate the displacement of the circle in the x, y
% coordinates from the origin

xcir = radius*cos(theta) + xdisp; % x coordinates of circle
ycir = radius*sin(theta) + ydisp; % y coordinates of circle

p2 = plot(xcir, ycir, 'b', 'linewidth', 1.8); % plot circle

% Assign value to do
do = radius;
disp(sprintf('The radius, thus the initial value of depth of magma chamber, of the encircled study area is do = %d m', do))
do = -5000;

%%%% Initial value of Xs and Ys %%%%

% Find weighted for x and y coordinates of Subsidence below a threshold

xWmean = abs(deform(subsid))'*procrd(subsid,1)/sum(abs(deform(subsid)));
yWmean = abs(deform(subsid))'*procrd(subsid,2)/sum(abs(deform(subsid)));

p3 = plot(xWmean, yWmean, '.r', 'markersize', 30); % plot the point
c2 = colorbar;
c2.Label.String = 'm/y';
c2.Label.FontSize = 13;

legend([p1, p2, p3], {'First assumption', 'Second assumption', 'Weighted mean'}, 'fontsize', 15, 'location', 'northwest')


%%%% Find initial value of DV %%%%
u = 0.27; % Poisson's ratio

DVn = (deform.*pi*do^2./(1 - u)).*(1 + ((procrd(:, 1) - xWmean).^2 + (procrd(:, 2) - yWmean).^2)./do^2).^(3/2);

DVo = mean(DVn);

clear minx maxx miny maxy xcirc ycirc
%%

%%%%%%%%%% Question 4 %%%%%%%%%%

%%%%% Symbolical computation of equations %%%%%

% Define symbolic variables
% p = pi, w = x, z = y, uP = poisson's ratio, xs = xso, ys = yso, d = do,
% DV = DVo
syms p w z uP xs ys d DV

% Define deformation (u) symbolic expression
def = ((1 - uP)*DV/(p*d^2))/(1 + ((w - xs)^2 + (z - ys)^2)/d^2)^(3/2);

k = [DV, d, xs, ys]; % row vector of unknown parameters
j = jacobian(def, k); % calculate jacobian w.r.t. uknow

j = matlabFunction(j); % convert symbolic expression j to matlabFunction

% Convert symbolic u equation to matlabFunction
def = matlabFunction(def);

%%

clear Du term1 term2 Dx

%%%%%%%%%% Question 5 %%%%%%%%%%

%%%%% Gauss-Newton implementation %%%%%

% Values initialization

% Initial values I chose at the previous Questions
ini = [DVo; do; xWmean; yWmean];
%ini = [ini zeros(4,100)];
Dx = ones(4,1); % increment step for every unknown
ini = [ini zeros(4,100)];
% Compose the iteration process

for i=1:100
    
    % Evaluate deformation function for current values of the unknowns
    uo = feval(def, ini(1, i), ini(2, i), pi, u, procrd(:,1), ini(3, i), ini(4,i), procrd(:,2));    
    
    % calculate observations minus evaluated deformation function for
    % current values for the unknown parameters
    Du = deform - uo;
    
    % Evaluate jacobian matrix for current values of the unknowns
    jacob = feval(j, ini(1, i), ini(2, i), pi, u, procrd(:,1), ini(3, i), ini(4, i), procrd(:,2));
    
    term1 = jacob'*Du;
    term2 = jacob'*jacob;
    
    Dx = term2\term1;
    ini(:, i+1) = ini(:, i) + Dx;
    
    % Condition to stop the iteration procedure
    if norm(Dx)^2 <= 1e-10
        break
    end
end

% clear the matrix with the estimated values from the zeros, based on the
% condition that every row of one column is equal to zero
ini(:, ini(1,:) == 0) = [ ];


%%%%% Plot unknowns and number of iteration %%%%%
[~, col] = size(ini);
figure
subplot(2,2,1)
plot(1:col, ini(1, :), '.r', 'markersize', 15);
xlabel('Number of iteration', 'fontsize', 13)
ylabel('Volume change (m^3)', 'fontsize', 13)
grid on

subplot(2,2,2)
plot(1:col, ini(2, :), '.b', 'markersize', 15);
xlabel('Number of iteration', 'fontsize', 13)
ylabel('Depth (m)', 'fontsize', 13)
grid on

subplot(2,2,3)
plot(1:col, ini(3, :), '.g', 'markersize', 15);
xlabel('Number of iteration', 'fontsize', 13)
ylabel('x coordinate (m)', 'fontsize', 13)
grid on

subplot(2,2,4)
plot(1:col, ini(4, :), '.m', 'markersize', 15);
xlabel('Number of iteration', 'fontsize', 13)
ylabel('y coordinate (m)', 'fontsize', 13)
grid on

clear row

%%

%%%%% Question 5 (Continue) %%%%%

figure
c = scatter(procrd(:,1), procrd(:,2), 15, deform, 'filled');
xlabel('X coordinates (m)', 'fontsize', 13);
ylabel('Y coordinates (m)', 'fontsize', 13);
title('Study area', 'fontsize', 16);
c = colorbar;
c.Label.String = 'm/y';
c.Label.FontSize = 13;
colormap(flipud(parula));
grid on
axis square
hold on

% xs = ini(3,34) and ys = ini(4,34) indicate the displacement of the circle in the x, y
% coordinates from the origin, based on the final estimates

xcir = radius*cos(theta) + ini(3,34); % x coordinates of circle
ycir = radius*sin(theta) + ini(4, 34); % y coordinates of circle

plot(xcir, ycir, '--r', 'linewidth', 4); % plot circle

plot(ini(3, 34), ini(4, 34), '.r', 'markersize', 25)

%%

%%%%%%%%%% Question 6 %%%%%%%%%%

sigma = 0.2e-3; % m/y converted to meters per year
J = feval(j, ini(1, 34), ini(2, 34), pi, u, procrd(:,1), ini(3, 34), ini(4, 34), procrd(:,2));

Q = sigma^2*inv(J'*J); % (co)variance matrix of all parameters
stds = sqrt(diag(Q)); % standard devations of every estimated parameter

CV = abs(stds*100./ini(:,34));
%%

%%%%%%%%%% Question 7 %%%%%%%%%%

% Calculate deformation rates for every x and y using the estimated DV, d, xs, ys

deformEst = feval(def, ini(1, 34), ini(2, 34), pi, u, procrd(:,1), ini(3, 34), ini(4, 34), procrd(:,2));

% Calculate residuals

resid = deform - deformEst;

figure
% Scatter of xy versus deformation
s1 = subplot(2,2,1);
scatter(procrd(:,1), procrd(:,2), 13, deformEst, 'filled');
xlabel('X coordinates (m)', 'fontsize', 13)
ylabel('Y coordinates (m)', 'fontsize', 13)
title('Modelled deformation spatial pattern', 'fontsize', 15);
c1 = colorbar;
c1.Label.String = 'm/y';
c1.Label.FontSize = 13;
colormap(s1, flipud(parula));
grid on

% Modelled deformation empirical PDF superimposing normal PDF
subplot(2,2,3)
histfit(deformEst, 20)
xlabel('Modelled deformation (m/y)', 'fontsize', 13)
ylabel('Density', 'fontsize', 13)
title('Empirical PDF of modelled Deformation', 'fontsize', 15);
grid on

% Scatter of xy versus residuals
s2 = subplot(2,2,2);
scatter(procrd(:,1), procrd(:,2), 13, resid, 'filled');
xlabel('X coordinates (m)', 'fontsize', 13)
ylabel('Y coordinates (m)', 'fontsize', 13)
title('Residuals spatial pattern', 'fontsize', 15);
c2 = colorbar;
c2.Label.String = 'm/y';
c2.Label.FontSize = 13;
colormap(s2, flipud(jet));
grid on

% Residual empirical PDF superimposing normal PDF
subplot(2,2,4)
histfit(resid, 20)
xlabel('Residuals (m/y)', 'fontsize', 13)
ylabel('Density', 'fontsize', 13)
title('Empirical PDF of Residuals', 'fontsize', 15);
grid on

%%

%%%%%%%%%%% Question 9 %%%%%%%%%%

% Calculate (co)variance matrix for every observation based on the
% (co)variance matrix generated for the estimates (DV, d, xs, ys)

Qu = J*Q*J';
ustds = sqrt(diag(Qu)); % calculate standard deviation for each point

% plot stds
figure
scatter(procrd(:,1), procrd(:,2), 15, ustds, 'filled')
xlabel('X coordinates (m)', 'fontsize', 13);
ylabel('Y coordinates (m)', 'fontsize', 13);
title('Deformation ó spatial pattern', 'fontsize', 15);
c2 = colorbar;
c2.Label.String = 'm/y';
c2.Label.FontSize = 13;
grid on
axis square