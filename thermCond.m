% Liquid Thermal Conductivity Function
%       Curvefit from the global variable 'F.k_fit' is needed.
% Ayaaz Yasin - Sep 11, 2024
%       k = thermCond(T)

function k = thermCond(T)
    global F
    k = polyval(F.k_fit, T);              % calculate k at given temperature(s)
end