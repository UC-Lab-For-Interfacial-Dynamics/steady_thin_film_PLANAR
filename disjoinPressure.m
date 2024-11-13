% Disjoining Pressure Function
% Reads Hamaker constant from 'fluidData.mat' and returns disjoining 
% pressure as a function of the input film height.
% Ayaaz Yasin - Sep 11, 2024
% Pd = disjoinPressure(h)

function Pd = disjoinPressure(h)
global F
% load('fluidData.mat', 'A')      % loads hamaker constant
Pd = F.A/(h^3);                   % disjoining pressure
end