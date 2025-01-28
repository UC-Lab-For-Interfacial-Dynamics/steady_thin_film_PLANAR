% Curvature Function
%       Calculates surface curvature as a function of first and second 
%       derivatives of the film height.
% Ayaaz Yasin - Sep 11, 2024
%
%       kappa = curvature(h1,h2)

function kappa = curvature(h1,h2)
    kappa = h2/(1 + h1^2)^(3/2);    % curvature in 2D planar coordinate system 
end