% Young Laplace Function
% Composite function to calculate and return vapor, disjoining, capillary,
% and liquid pressures.
% Ayaaz Yasin - Sep 11, 2024
% [Pv, Pl, Pc, Pd] = youngLaplace(T, h, h1, h2, h3)

function [Pv, Pl, Pc, Pd] = youngLaplace(T, h, h1, h2, h3)
Pv = vapPressure(T);                % call vapor pressure function
Pd = disjoinPressure(h);            % call disjoining pressure function
Pc = capPressure(T, h, h1, h2, h3); % call capillary pressure function
Pl = Pv - Pc - Pd;                  % calculate liquid pressure using the Young-Laplace equation