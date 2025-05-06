% Wall temperature function
%       Polynomial fit of the wall tempeature from MD.
% Ayaaz Yasin - Apr 18, 2025
%
%       Tw = Tw_func_MD(x)
%
% Input:    x_in .................. location along the wall [meters]
% Output:   Tw .................... vapor temperature [K]

function Tw = Tw_func_MD_16(x)
    x = x*1e9 + 14.4;   % offsetting input and converting to nm
    coeff = [-2.18061027e-02, 1.25124187e+00, 1.33661153e+02 ];
    Tw = polyval(coeff, x);

end