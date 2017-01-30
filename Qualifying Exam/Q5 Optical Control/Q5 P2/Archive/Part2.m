% Tim Coon - 9 December, 2014
% Qualifying Exam, Question #5, Part 2
% Introduce the control actuators to minimize the focal point distribution
% as a result of the random disturbance. Because the distribution is
% random, the evolution of the control is a random variable. What is the
% nature of the distribution of this control random variable? In
% particular, is there any bias and, if so, how does the noise distribution
% influence that bias?
clear all; close all; clc;

global f

% generate sample realizations and determine statistics of resulting focal-
% point-evolution ball radius
nreals = 1;
radius = zeros(nreals,1);
for i = 1:nreals
    [t,z,radius(i)] = Q5_Realization();
end

meanRadius = mean(radius);
stdRadius = std(radius);

%% Plot the last realization

x      = z(:,1);
xd     = z(:,2);
y      = z(:,3);
yd     = z(:,4);
alpha  = z(:,5);
alphad = z(:,6);

% focal point position
xfoc_nom = f*cos(pi/4);
yfoc_nom = f*sin(pi/4);
xfoc = x+f*cos(alpha+pi/4);
yfoc = y+f*sin(alpha+pi/4);

figure()
% suptitle('Motion Plots')
subplot(221)
plot(x,y)
title('Vertex Motion')
axis equal
subplot(222)
plot(xfoc,yfoc)
hold on
tc = linspace(0,2*pi,50);
xc_m = meanRadius*cos(tc)+xfoc_nom;
yc_m = meanRadius*sin(tc)+yfoc_nom;
plot(xc_m,yc_m)
xc_std = stdRadius*cos(tc)+xfoc_nom;
yc_std = stdRadius*sin(tc)+yfoc_nom;
plot(xc_std,yc_std,'--k')
hold off
xlim([xfoc_nom-2*meanRadius, xfoc_nom+2*meanRadius]);
ylim([yfoc_nom-2*meanRadius, yfoc_nom+2*meanRadius]);
axis equal
title('Focal Point Motion')
subplot(2,2,3:4)
plot(t,alpha)
title('Mirror Angle Change')

% arrow setup
Start = [x,y];
Stop  = [xfoc,yfoc];

% figure()
% arrow(Start,Stop,'Length',0.3,'BaseAngle',15,'TipAngle',25)
% axis square; axis equal

