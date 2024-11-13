% main - steady thin film setup file - PLANAR coordinate system
% Ayaaz Yasin - Nov 13, 2024
clear; clear global; clc; close all;
load fluidData.mat
global F C R    % global data sharing. F-fluid properties, C-constants, R-report/solution data

%%% Constants
[C.a,~] = tst_alpha(polyval(F.rhov_fit,140),polyval(F.rhol_fit,140));  % [-],          accommodation coefficient
C.Rc      = 8.31446261815324;         % [J/mol-K],    universal gas constant
C.Tv      = 140;                     % [K],          vapor temperature
Tw1       = 145;                    % [K],          wall temperature at bulk
Tw2       = 145;                    % [K],          wall temperature at adsorbed end
C.Tw_func = @(x) (Tw2-Tw1)/150e-6*x + Tw1; % [K],          wall temperature function

%%% Boundary conditions
% Ti      = C.Tw;                     % [K],          initial interface temperature
h         = 4e-9;                      % [m],          film height
h1        = -10;                      % [-],          first derivative
h2        = 1e2;                        % [m^-1],       second derivative
IC        = [h, h1, h2];                % [-],          initial conditions vector
C.x_in    = 0;                        % [m],          initial wall distance
C.L       = 50e-9;                  % [m],          length of solution
C.dx      = 1e-12;                     % [m],          step size

Ti = precondition(h,h1,h2);           % [K],          initial interface temperature

%%% Initialize report vectors
R.n        = 1;
R.xEval    = C.x_in;
R.h3       = nan;
R.mflux    = nan;
R.Ti       = Ti;
R.Tw       = C.Tw_func(C.x_in);
R.kappa    = curvature(h1,h2);
R.Pc       = nan;
R.Pd       = nan;
R.cv1 = nan;
R.cv2 = nan;

%% ODE solution
fprintf('distance fraction\t|\tfilm height\n------------------------------------\n')
[X, H] = ode4(@(X, H, step) evolution(X, H, step), [C.x_in, C.x_in+C.L], IC', C.dx);
H = H';

%%% Post-Processing
%% plots
close all;
fig=figure; fig.Position=[-1888,179,1815,706];
rmv = 10;

subplot(2,4,1); % film height
plot(X(1:end-rmv)*1e9, H(1:end-rmv,1)*1e9, '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [nm]'); ylabel('film height (h) [nm]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.YAxis(1).Exponent=0; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,2); % first derivative
plot(X(1:end-rmv)*1e9, H(1:end-rmv,2), '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [nm]'); ylabel('first derivative (h1) [-]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,3); % second derivative
plot(X(1:end-rmv)*1e9, H(1:end-rmv,3), '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [nm]'); ylabel('second derivative (h2) [1/m]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,4); % third derviative
plot(X(1:end-rmv)*1e9, R.h3(1:end-rmv), '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [nm]'); ylabel('third derivative (h3) [1/m^2]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,5); % mass flux
plot(X(1:end-rmv)*1e9, R.mflux(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [nm]'); ylabel('mass flux [kg/m^3-s]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,6); % interface temperature
plot(X(1:end-rmv)*1e9, R.Ti(1:end-rmv),'r.-', 'LineWidth',1.5, 'MarkerSize',10); hold on
plot(X(1:end-rmv)*1e9, R.Tw(1:end-rmv),'k.-', 'LineWidth',1.5, 'MarkerSize',10);
legend('interface','wall','Location','best')
xlabel('length along wall (x) [nm]'); ylabel('temp [K]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,7); % curvature
plot(X(1:end-rmv)*1e9, R.kappa(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [nm]'); ylabel('curvature [1/m]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(2,4,8); % Pressures
plot(X(1:end-rmv)*1e9, R.Pc(1:end-rmv),'r.-', 'LineWidth',1.5, 'MarkerSize',10); hold on %ylabel('capillary pressure [Pa]');
plot(X(1:end-rmv)*1e9, R.Pd(1:end-rmv),'k.-', 'LineWidth',1.5, 'MarkerSize',10); %ylabel('disjoining pressure [Pa]');
xlabel('length along wall (x) [nm]'); ylabel('pressure [Pa]')
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded; 
legend('capillary','disjoining','location','best')