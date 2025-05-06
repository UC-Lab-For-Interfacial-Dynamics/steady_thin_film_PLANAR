% Main Setup and Post-Processing Script
%       Steady-state thin film solution in a 2D planar coordinate system.
%       See 'thin films notes.pdf' for details on the models.
%       This script works in three sections: 
%       (1) Numerical & Problem Setup - adjust initial/boundary conditions
%           here. Global variables are used to pass parameters to other
%           functions. Wall and vapor temps can be defined as function handles.
%       (2) ODE solution - the problem is passed to a Runge-Kutta solver.
%           If needed, change the solver function here. 
%       (3) Post-processing - mass flow integration and solution plots are generated.
%
% Ayaaz Yasin - Nov 13, 2024
%

clear; clear global; clc; close all;
warning('off','all')
load fluidData.mat      % load fluid properties in structure 'F' - see Ar_data_file_generator.m
global F C R            % global data sharing. F-fluid properties, C-constants, R-report/solution data

%%% Constants, Setup, and Input Parameters
C.debug   = 'off';                  % when 'on', plots inner loop iterations & prints ODE solution at every step.
C.Rc      = 8.31446261815324;       % [J/mol-K],    universal gas constant
C.Lz      = 3.92e-9;                % depth of MD domain [m]
C.Tv_func = @Tv_func_MD_centerline_16; % vapor temperature function handle
C.Tw_func = @Tw_func_MD_16;            % wall temperature function handle

%%% Boundary Conditions
C.x_in    = 0.1e-9;                 % [m],          initial wall distance
h         = film_16(C.x_in);        % [m],          film height (currently from MD)
h1        = -6.265107527173289e+00; % [-],          first derivative (currently from MD)
h2        = 3.192157326161009*1e9;  % [m^-1],       second derivative (currently from MD)
R.G       = 2.82584020795487e-14;   % [kg/s],       inlet mass flow rate
IC        = [h,h1,h2];              % [-],          initial conditions vector
C.L       = 30e-9;                  % [m],          length of solution
C.dx      = 1e-12;                  % [m],          step size

%%% Initialize report vectors
R.n        = 1;                     % iteration count
R.xEval    = C.x_in;                % x-vector
R.h3       = nan;                   % third derivative vector
R.Tw       = C.Tw_func(C.x_in);     % wall temp vector
R.Tv       = C.Tv_func(C.x_in);     % vapor temp vector
[R.Ti, R.mflux] = precondition(h,h1,h2, R.Tv, R.Tw);   % [K], initial interface temperature
R.kappa    = curvature(h1,h2);      % curvature vector
R.Pc       = surfTension(R.Ti)*R.kappa;   % capillary pressure vector
R.Pd       = disjoinPressure(h);    % disjoining pressure vector
R.a        = nan;                   % accommodation coefficient vector
R.mflow    = 0;                     % evaporation mass flow vector
R.pdGrad   = nan;                   % disjoining pressure gradient vector
R.plGrad   = nan;                   % liquid pressure gradient vector
R.stGrad   = nan;                   % surface tension gradient vector
R.h        = h;                     % film height vector
R.surfArea = 0;                     % cummulative surface area
R.totMflow = 0;                     % cummulative phase change mass flow


%% ODE solution
[X, H] = ode4(@(X, H, step) evolution(X, H, step), [C.x_in, C.L], IC', C.dx);
H = H'; % transposing for standard indexing

%% Post-Processing
%%% converting mass flux to mass flow

mflow = massIntegration(X, H(:,1), R.mflux);
fprintf('Total phase change mass flux %0.4e kg/s\t\t%0.4e kg/s\n', mflow, sum(R.mflow))

%%% plots
close all;
fig=figure; fig.Position=[3,50,1917,946];
rmv = 2;        % number of points to remove from the end of the solution in the plots. Useful if the solution blows up.

subplot(3,4,1); % film height
plot(X(1:end-rmv)*1e10, H(1:end-rmv,1)*1e10, 'r.-', 'LineWidth',1.5, 'markersize',10); hold on;
plot(X(1:end-rmv)*1e10, film_16(X(1:end-rmv))*1e10, 'k--', 'LineWidth',1.5);
xlabel('length along wall (x) [A]'); ylabel('film height (h) [A]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.YAxis(1).Exponent=0; ax.XAxis(1).Exponent=0; axis padded;
legend('ODE solution','MD','location','best')

subplot(3,4,2); % first derivative
plot(X(1:end-rmv)*1e10, H(1:end-rmv,2), '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [A]'); ylabel('first derivative (h1) [-]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,3); % second derivative
plot(X(1:end-rmv)*1e10, H(1:end-rmv,3), '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [A]'); ylabel('second deriv (h2) [1/m]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,4); % third derviative
plot(X(1:end-rmv)*1e10, R.h3(1:end-rmv), '.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [A]'); ylabel('third derivative (h3) [1/m^2]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,5); % mass flux & accommodation coefficient
yyaxis left 
    plot(X(1:end-rmv)*1e10, R.mflux(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10);
    ylabel('mass flux [kg/m^2-s]');
yyaxis right
    plot(X(1:end-rmv)*1e10, R.a(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10);
    ylabel('accommodation coeff [-]');
xlabel('length along wall (x) [A]'); 
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,6); % interface temperature
plot(X(1:end-rmv)*1e10, R.Ti(1:end-rmv),'r.-', 'LineWidth',1.5, 'MarkerSize',10); hold on
plot(X(1:end-rmv)*1e10, Ti_func_MD_16(X(1:end-rmv)),'r--', 'LineWidth',1);
plot(X(1:end-rmv)*1e10, R.Tw(1:end-rmv),'k.-', 'LineWidth',1.5, 'MarkerSize',10);
plot(X(1:end-rmv)*1e10, R.Tv(1:end-rmv),'b.-', 'LineWidth',1.5, 'MarkerSize',10);
legend('interface','interface MD','wall','vapor','Location','best')
xlabel('length along wall (x) [A]'); ylabel('temp [K]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,7); % curvature
plot(X(1:end-rmv)*1e10, R.kappa(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [A]'); ylabel('curvature [1/m]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,8); % Pressures
plot(X(1:end-rmv)*1e10, R.Pc(1:end-rmv),'r.-', 'LineWidth',1.5, 'MarkerSize',10); hold on 
plot(X(1:end-rmv)*1e10, R.Pd(1:end-rmv),'k.-', 'LineWidth',1.5, 'MarkerSize',10);
xlabel('length along wall (x) [A]'); ylabel('pressure [Pa]')
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded; 
legend('capillary','disjoining','location','best')

subplot(3,4,9); % mass flow plots
yyaxis left 
    plot(X(2:end-rmv)*1e10, R.mflow(2:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10);
    ylabel('phase change mass flow [kg/s]'); 
    ylim([min(R.mflow(2:end-rmv)),max(R.mflow(2:end-rmv))]);
yyaxis right
    plot(X(1:end-rmv)*1e10, R.G(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10);
    ylabel('liquid mass flow [kg/s]');
xlabel('length along wall (x) [A]'); 
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax.XAxis(1).Exponent=0; %axis padded;

subplot(3,4,10); % liquid pressure gradient
plot(X(1:end-rmv)*1e10, R.plGrad(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [A]'); ylabel('liquid pressure grad [Pa/m]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;

subplot(3,4,11); % disjoining pressure gradient
plot(X(1:end-rmv)*1e10, R.pdGrad(1:end-rmv),'.-', 'LineWidth',1.5, 'MarkerSize',10, 'MarkerEdgeColor', 'r');
xlabel('length along wall (x) [A]'); ylabel('disjoining pressure grad [Pa/m]');
set(gca, 'FontWeight', 'bold', 'fontsize', 14); ax = gca; ax.XAxis(1).Exponent=0; axis padded;
