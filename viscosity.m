% Viscosity Function
% Reads polyfit coefficients from 'fluidData.mat' and returns dynamic and 
% kinematic viscosities as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% [mu, nu] = viscosity(T)

function mu = viscosity(T)
global F
% load('fluidData.mat', 'mu_fit')     % loads polyfit coefficients for mu
mu = polyval(F.mu_fit, T);            % calculate mu at given temperature(s)
% rho = density(T);                   % call density function
% nu = mu./rho;                       % kinematic viscosity
end