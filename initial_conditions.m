clear; clc; close all;

x = 245.6; %[A]
h = film(0);
a = 210.6034273381962;
b = 37.43510914754422;
c = 35;

h1 = -(x-a)/(h-b);
h2 = ((x-a)^2 - (h-b)^2)/((h-b)^3);

fprintf('%0.16e\n%0.16e\n', h1, h2)