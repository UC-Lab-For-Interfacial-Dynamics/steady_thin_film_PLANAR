% Liquid Thermal Conductivity Function
% Reads polyfit coefficients from 'fluidData.mat' and returns liquid
% thermal conductivity as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% k = thermCond(T)

function k = thermCond(T)
global F
% load('fluidData.mat', 'k_fit')      % loads polyfit coefficients for thermal conductivity
k = polyval(F.k_fit, T);              % calculate k at given temperature(s)
end