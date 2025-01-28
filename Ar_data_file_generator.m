% Argon data file generator
% Script file that reads NIST data from Excel files, generates curvefits,
%       and generates fluidData.mat containing a structure 'F' that can be 
%       read by the numerical solution code. Modify this script to change
%       fluid properties or the fluid.
% Ayaaz Yasin - Nov 13, 2024
clear; clc; close all;

data = importdata('NIST_Ar_saturated_liquid_100K-150K.xlsx');
data = data.data;
if isempty(data); error('data is empty. excel file might be open.'); end

T       = data(:,1);        % [K]       saturation temperature
P       = data(:,2)*1e6;    % [Pa]      saturation pressure
rhol    = data(:,3);        % [kg/m^3]  liquid density
V       = data(:,4);        % [m^3/kg]  volume
h_liq   = data(:,6)*1e3;    % [J/kg]    liquid enthalpy
cp      = data(:,9)*1e3;    % [J/kg-K]  specific heat
mu      = data(:,12);       % [Pa-s]    dynamic viscosity
k       = data(:,13);       % [W/m-K]   thermal conductivity 
sigma   = data(:,14);       % [N/m]     surface tension

% vapor properties
data = importdata('NIST_Ar_saturated_vapor_100K-150K.xlsx');
data = data.data;
rhov        = data(:,3);        % [kg/m^3]  density
h_vap       = data(:,6)*1e3;    % [J/kg]    vapor enthalpy

% calculated properties
nu      = mu./rhol;          % [m^2/s]   kinetmatic viscosity
hfg     = h_vap - h_liq;     % [J/kg]    enthalpy of vaporization

%% generate polyfits
F.Psat_fit    = polyfit(T, P,     3);
F.rhol_fit    = polyfit(T, rhol,  5);
F.rhov_fit    = polyfit(T, rhov,  5);
F.V_fit       = polyfit(T, V,     6);
F.cp_fit      = polyfit(T(1:end-3), cp(1:end-3),    6);
F.mu_fit      = polyfit(T, mu,    4);
F.k_fit       = polyfit(T, k,     3);
F.sigma_fit   = polyfit(T, sigma, 1);
F.hfg_fit     = polyfit(T, hfg,   5);

%% setup fluidData.mat
F.name    = 'argon';
F.M       = 0.039948;       % [kg/mol]      molar mass
F.A       = 9.91e-21;       % [J]           hamaker/dispersion constant (from MD data)

save('fluidData.mat', 'F');

