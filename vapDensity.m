% Vapor Density Function
% Reads polyfit coefficients from 'fluidData.mat' and returns vapor
% density as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% rhov = vapDensity(T)

function rhov = vapDensity(T)
load('fluidData.mat', 'rhov_fit')     % loads polyfit coefficients for density
rhov = polyval(rhov_fit, T);          % calculate rhov at given temperature(s)
end