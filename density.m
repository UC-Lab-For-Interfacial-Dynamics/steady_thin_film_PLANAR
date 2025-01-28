% Liquid Density Function
%       Curvefit from the global variable 'F.rhov_fit' and 'F.rhol_fit' are needed.
% Ayaaz Yasin - Sep 11, 2024
%
%       [rhov, rhol] = density(T)
%
% Input:    T ............ Temperature [K]
% Outputs:  rhov ......... vapor density [kg/m^3]
%           thol ......... liquid density [kg/m^3]

function [rhov, rhol] = density(T)
    global F
    rhov = polyval(F.rhov_fit, T);   % calculate vapor rho at given temperature(s)
    rhol = polyval(F.rhol_fit, T);   % calculate liquid rho at given temperature(s)
end