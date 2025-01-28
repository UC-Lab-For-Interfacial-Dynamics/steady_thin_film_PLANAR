% Vapor Saturation Pressure Function
%       Curvefit from the global variable 'F.Psat_fit' is needed.
% Ayaaz Yasin - Sep 11, 2024
% Pv = vapPressure(T)

function Pv = vapPressure(T)
    global F
    Pv = polyval(F.Psat_fit, T);           % calculate Pv at given temperature(s)
end