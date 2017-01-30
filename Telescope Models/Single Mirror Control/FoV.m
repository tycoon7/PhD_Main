% plot FoV vs focal length
% http://en.wikipedia.org/wiki/Angle_of_view
clear all; close all; clc;

f = 0.7:0.001:0.8;
d = 0.5;

alpha = 2*atand(d/2*f);

figure(1)
plot(f,alpha)
title('FoV Characterization')
xlabel('Focal Length (m)'); ylabel('Half-FoV (deg)')