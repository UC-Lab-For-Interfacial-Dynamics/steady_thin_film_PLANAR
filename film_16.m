% SMU MD film height piecewise curvefits for 16 nm
% converted to
%   (1) all dimensions in meters
%   (2) apex of bulk at x=0, thin and adsorbed films are towards positive x.
% Ayaaz Yasin - Apr 18, 2025
%
%       h = film(X)
%
% Input:   x ............ along channel length [m]
% Output:  h ............ film height [m]

% translational offset between UC and SMU coordinate systems
% x = 0 in UC frame is 14.4 nm in SMU's frame

function H = film_16(X)
X = X*1e9 + 14.4;   % offsetting input and converting to nm
H = [];
for x = X
if 14.4<=x && x<=14.4985
    a = 22.4;
    b = 8.428;
    r = 8;
    h = b - sqrt(r^2 - (x-a)^2);
elseif 14.4985<x && x<=15.8
    coeff = [2.4393724485, -149.6432063298, 3442.2657550721, -35193.1369557010, 134944.6860045021];
    h = polyval(coeff, x);
elseif 15.8<x && x<=20
    coeff = [-0.0021121337, 0.2060672634,-8.0299573876,156.2357093721,-1517.9881051692,5894.6497277090];
    h = polyval(coeff, x);
elseif 20<x && x<=31.2
    coeff = [-0.0130684023,1.7098974910];
    h = polyval(coeff, x);
else
    h = nan;
    % error('invalid x input. no film exists here.')
end

h = h + 0.428; % wall offset
h = h*1e-9; % convert to meters
H(end+1) = h;
end
end