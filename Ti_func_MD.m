% Interface temperature function
%       Polynomial fit of the interface tempeature from MD.
% Ayaaz Yasin - Jan 30, 2025
%
%       Ti = Ti_func_MD(x_in)
%
% Input:    x_in .................. location along the interface [meters]
% Output:   Ti .................... vapor temperature [K]

function Ti = Ti_func_MD(x_in)

% polynomial fit from SMU. change this if needed.
Ti_MD = @(x) -5.19492236e-07*x.^4 + 3.68657055e-04*x.^3 -1.00873008e-01*x.^2 +1.25653037e+01*x  -4.55367990e+02;

% transformation to meters and flipping the x-axis so the meniscus apex is at x=0 and reducing film film height at increasing x.
Ti_func = @(x) Ti_MD(-(1e10*x-245.6));

Ti = Ti_func(x_in);

end