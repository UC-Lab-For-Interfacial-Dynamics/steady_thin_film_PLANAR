% Specfic Heat Function
% Reads polyfit coefficients from 'fluidData.mat' and returns liquid
% specific heat as a function of the input temperature.
% Ayaaz Yasin - Sep 11, 2024
% cp = specificHeat(T, cp_fit)

function cp = specificHeat(T, cp_fit)
% load('fluidData.mat', 'cp_fit')    % loads polyfit coefficients for specific heat
cp = polyval(cp_fit, T);           % calculate cp at given temperature(s)
end