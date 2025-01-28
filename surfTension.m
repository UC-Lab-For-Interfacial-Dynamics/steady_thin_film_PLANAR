% Liquid Density Function
%       Value of surface tension is currently set to a negative constant from 
%       the MD data. Uses curvefit from NIST data to find the surface tension 
%       gradient with temperature (ds/dT).
% Ayaaz Yasin - Sep 11, 2024
%       [sigma, gamma] = surfTension(T)

function [sigma, gamma] = surfTension(T)
    global F
    % sigma = polyval(F.sigma_fit, T);   % uncomment to calculate sigma using NIST curvefit at given temperature(s)
    sigma = -3.7e-3;                   % [N/m], surf tension from MD data
    gamma = F.sigma_fit(1);            % surface tension gradient
end