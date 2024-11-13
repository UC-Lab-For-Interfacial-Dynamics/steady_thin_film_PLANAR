% Specfic Volume Function
% Reads polyfit coefficients from 'fluidData.mat' and returns liquid
% specific volume as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% Vl = specificVol(T)

function Vl = specificVol(T)
global F
% load('fluidData.mat', 'V_fit')    % loads polyfit coefficients for specific volume
Vl = polyval(F.V_fit, T);           % calculate Vl at given temperature(s)
end