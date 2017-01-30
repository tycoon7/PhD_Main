% parametric circle and parabola

clear; close all; clc;

t = linspace(0,1);
f = 1;

x_p = -2*f*t;
y_p = f*t.^2;

x_c = -(1-t.^2)./(1+t.^2);
y_c = 1 - (2*t./(1+t.^2));

figure()
plot(x_p,y_p,x_c,y_c);
% plot(x_c,y_c);
axis equal
