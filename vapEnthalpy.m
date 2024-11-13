% Enthalpy of Vaporization Function
% Reads polyfit coefficients from 'fluidData.mat' and returns hfg as a 
% function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024

function hfg = vapEnthalpy(T)
load('fluidData.mat', 'hfg_fit')    % loads polyfit coefficients for hfg
hfg = polyval(hfg_fit, T);          % calculate hfg at given temperature(s)
end