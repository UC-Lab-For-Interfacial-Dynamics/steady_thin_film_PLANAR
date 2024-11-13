% Curvature Function
% Calculates surface curvature and curvature derivature as a function of
% first, second, and third derivatives of the film height.
% Ayaaz Yasin - Sep 11, 2024
% kappa = curvature(h,h1,h2)

function kappa = curvature(h1,h2)
global C

% curvature in the 2D planar coordinate system
term1 = h2/(1 + h1^2)^(3/2);
term2 = 0;
kappa =  term1+term2;

end