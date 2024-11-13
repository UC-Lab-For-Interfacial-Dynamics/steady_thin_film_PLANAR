% Surface Tension Function
% Reads polyfit coefficients from 'fluidData.mat' and returns surface
% tension and surface tension gradient as a function of the input 
% temperature.
% Ayaaz Yasin - Sep 11, 2024
% [sigma, gamma] = surfTension(T)

function [sigma, gamma] = surfTension(T)
global F
% load('fluidData.mat', 'sigma_fit')    % loads polyfit coefficients for surface tension
sigma = polyval(F.sigma_fit, T);        % calculate sigma at given temperature(s)
gamma = F.sigma_fit(1);                 % surface tension gradient
end