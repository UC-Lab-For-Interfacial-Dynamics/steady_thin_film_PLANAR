% mass flux integration function
%       function converts from phase change mass flux profile (kg/m^2-s) to 
%       the total mass flow rate (kg/s). Using the film profile.
% Ayaaz Yasin - Jan 27, 2025
%
%       mflow = massIntegration(x, h, mflux)
%
% Inputs:   (all vectors must be the same size)
%       x .................. vector of horizontal grid points [m]
%       h .................. film height vector [m]
%       mflux .............. phase change mass flux vector [kg/m^2-s]
% Output:
%       mflow .............. phase change mass flow [kg/s]

function    mflow = massIntegration(x, h, mflux)
global C
    mflow = 0;
    for i = 1:length(x)-1
        dA  = C.Lz*sqrt((h(i)-h(i+1))^2 + C.dx^2);    % area element  [m^2]
        mflow = mflow + dA*mflux(i);
    end
end