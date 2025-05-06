clear; clc; close all;

xstart = 0.1e-9;   % starting x is offset to avoid large derivative
x = linspace(xstart,(31.2-14.4)*1e-9, 1e6); 
H = nan(size(x));
for i = 1:length(x)
    H(i) = film_16(x(i));
end
figure; 
yyaxis left
plot(x,H,'.'); hold on; axis equal;


% dx = x(2)-x(1);
% h1 = ((-49/20)*H(1) + (6)*H(2) + (-15/2)*H(3) + (20/3)*H(4) + (-15/4)*H(5) + (6/5)*H(6) + (-1/6)*H(7))/dx;
% h2 = ((469/90)*H(1) + (-223/10)*H(2) + (879/20)*H(3) + (-949/18)*H(4) + (41)*H(5) + (-201/10)*H(6) + (1019/180)*H(7) + (-7/10)*H(8))/dx^2;

x1 = xstart*1e9 + 14.4;
a = 22.4;
b = 8.428;
r = 8;
h = H(1)
h1= (x1-a)/sqrt(r^2 - (x1-a)^2)
h2= (r^2)/(r^2 - (x1-a)^2)^(3/2)

Tv = Tv_func_MD_centerline_16(x);
yyaxis right
plot(x, Tv, '.')