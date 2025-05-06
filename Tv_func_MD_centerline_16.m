% Vapor temperature function
%       Polynomial fit of the centerline vapor tempeature from MD.
% Ayaaz Yasin - Apr 18, 2025
%
%       Tv = Tv_func_MD_centerline(x_in)
%
% Input:    x_in .................. location along the wall [meters]
% Output:   Tv .................... vapor temperature [K]


function Tv = Tv_func_MD_centerline_16(x)
    x = x*1e9 + 14.4;   % offsetting input and converting to nm
    coeff = [-0.0003117, 0.03116, -1.2097, 22.55,-28.81];
    Tv = polyval(coeff, x);
end