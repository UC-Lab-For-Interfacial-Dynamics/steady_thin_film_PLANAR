% Capillary Pressure Function
% Calculates capillary pressure and using surface tension and surface
% curvature by calling their respective functions.
% Ayaaz Yasin - Sep 11, 2024
% Pc = capPressure(T, h, h1, h2, h3, sigma_fit)

function Pc = capPressure(T, h, h1, h2, h3)
global F C
[sigma, ~] = surfTension(T, F.sigma_fit);                % call surface tension function
[kappa, ~] = curvature(h1, h2);   % call surface curvature function
Pc = sigma*kappa;
end