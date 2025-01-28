% mass flux integration function
%       function converts from phase change mass flux (kg/m^2-s) to the 
%       mass flow rate (kg/s). Using the film profile.
% Ayaaz Yasin - Jan 27, 2025
%
%       mflow = massIntegration(x, h, mflux)
%
% Inputs:
%       x .................. vector of horizontal grid points [m]
%       h .................. film height vector [m]
%       mflux .............. phase change mass flux vector [kg/m^2-s]
%       all vectors must be the same size.
% Output:
%       mflow .............. phase change mass flow [kg/s]
%                            will be 1 element shorter than input vectors

function    mflow = massIntegration(x, h, mflux)

mflow = nan(size(x));
for i = 1:length(x)-1
    dx = x(i+1)-x(i);                       % in case of non-uniform solution [m^2]
    A  = sqrt((h(i)-h(i+1))^2 + dx^2);      % area element  [m^2]
    mflow(i) = A*mflux(i);
end

mflow = sum(mflow(1:end-1));
end