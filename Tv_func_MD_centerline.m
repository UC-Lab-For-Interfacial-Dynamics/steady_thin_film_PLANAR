% Vapor temperature function
%       Polynomial fit of the centerline vapor tempeature from MD.
% Ayaaz Yasin - Jan 26, 2025
%
%       Tv = Tv_func_MD_centerline(x_in)
%
% Input:    x_in .................. location along the wall [meters]
% Output:   Tv .................... vapor temperature [K]


function Tv = Tv_func_MD_centerline(x_in)

    % polynomial fit from SMU. change this if needed.
    Tv_MD = @(x) -0.00236*x.^2 + 0.73094*x + 83.23931;  
    
    % transformation to meters and flipping the x-axis so the meniscus apex is at x=0 and reducing film film height at increasing x.
    Tv_func = @(x) Tv_MD(-(1e10*x-245.6));  
    
    Tv = Tv_func(x_in);

end