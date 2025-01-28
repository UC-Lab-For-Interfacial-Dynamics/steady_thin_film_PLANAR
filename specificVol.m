% Specfic Volume Function
%       Curvefit from the global variable 'F.V_fit' is needed.
% Ayaaz Yasin - Sep 11, 2024
%       Vl = specificVol(T)

function Vl = specificVol(T)
    global F
    Vl = polyval(F.V_fit, T);           % calculate Vl at given temperature(s)
end