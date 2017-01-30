% Tim Coon - 9 December, 2014
% Qualifying Exam, Question #5, Part 2
% Introduce the control actuators to minimize the focal point distribution
% as a result of the random disturbance. Because the distribution is
% random, the evolution of the control is a random variable. What is the
% nature of the distribution of this control random variable? In
% particular, is there any bias and, if so, how does the noise distribution
% influence that bias?
clear; close all; clc;

global r0 sigma f

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

% focal point position
xfoc_nom = f*cos(pi/4);
yfoc_nom = f*sin(pi/4);
xfoc = x + f*cos(alpha+pi/4);
yfoc = y + f*sin(alpha+pi/4);

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
title('CoM Linear Disp')
subplot(2,2,4)
plot(t,alpha)
title('Angular Disp')


% arrow setup
% Start = [x,y];
% Stop  = [xfoc,yfoc];
% 
% figure()
% arrow(Start,Stop,'Length',0.3,'BaseAngle',15,'TipAngle',25)
% axis square; axis equal

%% Frequency analysis

dt = t(2);      % (s) sample time
Fs = 1/dt;      % (Hz) sample frequecy
Lx = length(x);
NFFT = 2^nextpow2(Lx); % Next power of 2 from length of x
X = fft(x,NFFT)/Lx;
freq = Fs/2*linspace(0,1,NFFT/2+1);

% Plot single-sided amplitude spectrum.
figure()
plot(freq,2*abs(X(1:NFFT/2+1))) 
title('Single-Sided Amplitude Spectrum of y(t)')
xlabel('Frequency (Hz)')
ylabel('|Y(f)|')
