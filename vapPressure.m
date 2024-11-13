% Vapor Saturation Pressure Function
% Reads polyfit coefficients from 'fluidData.mat' and returns saturation
% pressure as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% Pv = vapPressure(T)

function Pv = vapPressure(T)
global F
% load('fluidData.mat', 'Psat_fit')    % loads polyfit coefficients for saturation pressure
Pv = polyval(F.Psat_fit, T);           % calculate Pv at given temperature(s)
end