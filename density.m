% Liquid Density Function
% Reads polyfit coefficients from 'fluidData.mat' and returns liquid
% density as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% rho = density(T)

function rho = density(T)
global F
% load('fluidData.mat', 'rhol_fit')    % loads polyfit coefficients for density
rho = polyval(F.rhol_fit, T);          % calculate rho at given temperature(s)
end