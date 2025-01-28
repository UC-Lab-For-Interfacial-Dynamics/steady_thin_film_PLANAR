% Viscosity Function
%       Curvefit from the global variable 'F.mu_fit' is needed. Calculates
%       dynamic viscosity.
% Ayaaz Yasin - Sep 11, 2024
%       mu = viscosity(T)

function mu = viscosity(T)
    global F
    mu = polyval(F.mu_fit, T);   % [Pa-s] calculate dynamic viscosity (mu) at given temperature(s)
end