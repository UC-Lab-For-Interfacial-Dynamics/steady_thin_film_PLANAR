% SMU MD film height piecewise curvefits for 8 nm
% converted to
%   (1) all dimensions in meters
%   (2) apex of bulk at x=0, thin and adsorbed films are towards positive x.
% Ayaaz Yasin - Jan 20, 2024
%
%       h = film(X)
%
% Input:   x ............ along channel length [m]
% Output:  h ............ film height [m]


function h = film_8(X)
X = -(X.*1e10-245.6);     % converting to angstorm
h = [];
for x = X
    if 160<=x && x<=220
        % h(end+1) = 0.036735*x + 9.315535;
        h(end+1) = 0.036735*x + 0.691535;
    elseif 220<x && x<=242
        % h(end+1) = 7.14271267e-05*x.^4 - 0.0646471235*x.^3 + 21.9539327*x.^2 - 3315.12577*x + 187814.967;
        h(end+1) = 7.14271267e-05*x.^4 - 0.0646471235*x.^3 + 21.9539327*x.^2 - 3315.12577*x + 187806.343;
    elseif 242<x && x<=245.6
        % a = 210.6034273381962;
        % b = 46.05910914754422;
        % r = 35;
        a = 210.6034273381962;
        b = 37.43510914754422;
        r = 35;
        h(end+1) = b - sqrt(r^2-(x-a).^2);
    else
        h(end+1) = nan;
    end
end
h = h.*1e-10;
% yWall = 9e-10;
% h = h-yWall;
end