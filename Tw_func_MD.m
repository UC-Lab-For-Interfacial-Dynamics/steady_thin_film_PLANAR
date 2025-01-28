% Wall temperature function
%       Polynomial fit of the wall tempeature from MD.
% Ayaaz Yasin - Jan 26, 2025
%
%       Tw = Tw_func_MD(x_in)
%
% Input:    x_in .................. location along the wall [meters]
% Output:   Tw .................... vapor temperature [K]

function Tw = Tw_func_MD(x_in)

% polynomial fit from SMU. change this if needed.
Tw_MD = @(x) -8.605153181387093e-05*x.^3 + 0.047714442101985444*x.^2 - 8.802478773174819*x + 683.9205403131378;

% transformation to meters and flipping the x-axis so the meniscus apex is at x=0 and reducing film film height at increasing x.
Tw_func = @(x) Tw_MD(-(1e10*x-245.6));

Tw = Tw_func(x_in);

end