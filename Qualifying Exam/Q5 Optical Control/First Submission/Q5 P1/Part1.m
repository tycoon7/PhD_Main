% Tim Coon - 9 December, 2014
% Qualifying Exam, Question #5, Part 1
% With no active control, characterize the evolution of the focal point
% distribution in response to the random disturbance input.
clear; close all; clc;

global f

% Geometry
R = 0.5;        % (m) mirror radius
f = R/2;        % (m) mirror focal length
R1 = 2*R/pi;
R2 = R-2*R/pi;
rho = norm([R1,R2]);
rho1 = [-R2; R1];
rho2 = [R1; -R2];
phi1 = atan(R2/R1);
phi2 = atan(R1/R2);
cos_phi1 = cos(phi1);
sin_phi1 = sin(phi1);
cos_phi2 = cos(phi2);
sin_phi2 = sin(phi2);

% generate sample realizations and determine statistics of resulting focal-
% point-evolution ball radius
nreals = 10;
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

% mirror vertex position
x2_nom = rho2(1);
y2_nom = rho2(2);
x2  = x2_nom + x + rho*alpha*cos_phi2;
y2  = y2_nom + y + rho*alpha*sin_phi2;

% focal point position
xfoc_nom = x2_nom;
yfoc_nom = y2_nom + f;
xfoc = x2 + f*sin(alpha);
yfoc = y2 + f*cos(alpha);

figure()
% suptitle('Motion Plots')
subplot(221)
plot(x,y)
title('CoM Motion')
axis equal
subplot(222)
plot(xfoc,yfoc)
hold on
tc = linspace(0,2*pi,50);
xc_m = meanRadius*cos(tc)+xfoc_nom;
yc_m = meanRadius*sin(tc)+yfoc_nom;
plot(xc_m,yc_m)
xc_std1 = (meanRadius - stdRadius)*cos(tc) + xfoc_nom;
yc_std1 = (meanRadius - stdRadius)*sin(tc) + yfoc_nom;
xc_std2 = (meanRadius + stdRadius)*cos(tc) + xfoc_nom;
yc_std2 = (meanRadius + stdRadius)*sin(tc) + yfoc_nom;
plot(xc_std1,yc_std1,'--k',xc_std2,yc_std2,'--k')
hold off
xlim([xfoc_nom-2*meanRadius, xfoc_nom+2*meanRadius]);
ylim([yfoc_nom-2*meanRadius, yfoc_nom+2*meanRadius]);
axis equal
legend('Pos','Mean','Std')
title('Focal Point Motion')
subplot(223)
plot(t,x,t,y)
legend('x','y')
title('Linear Disp')
subplot(2,2,4)
plot(t,alpha)
title('Angular Disp')

% arrow setup
Start = [x,y];
Stop  = [xfoc,yfoc];

% figure()
% arrow(Start,Stop,'Length',0.3,'BaseAngle',15,'TipAngle',25)
% axis square; axis equal

