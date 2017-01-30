% Tim Coon - 9 December, 2014
% Qualifying Exam, Question #5, Part 1
% How does the distribution of the disturbance affect the distribution of
% the output?
clear; close all; clc;

global r0 sigma f

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
nreals = 100;
r0 = 0:0.1:0.5;
sigma = 0:0.02:0.1;
radius = zeros(nreals,1);
meanRadius = zeros(length(r0),length(sigma));
stdRadius = zeros(length(r0),length(sigma));
for i = 1:length(r0)
    for j = 1:length(sigma)
        for k = 1:nreals
            [t,z,radius(k)] = Q5_Realization(r0(i),sigma(j));       
        end
        meanRadius(i,j) = mean(radius);
        stdRadius(i,j) = std(radius);
    end
end

% save('OutputData_100_Realizations.mat');

%% Plot stats out wrt stats in

load('OutputData_100_Realizations.mat');

figure()
suptitle('Output Stats Vs. Input Stats')
subplot(221)
hold on
for i = 1:length(sigma)
    plot(r0,meanRadius(:,i),'DisplayName',['\sigma = ',num2str(sigma(i))]);
end
hold off
xlabel('Mean Disturbance Strength In'); ylabel('Mean Radius Out');
legend(gca,'show','Location','SouthEast')

% figure()
subplot(222)
hold on
for i = 1:length(sigma)
    plot(r0,stdRadius(:,i),'DisplayName',['r_0 = ',num2str(r0(i))])
end
hold off
xlabel('Mean Disturbance Strength In'); ylabel('\sigma Radius Out');
legend(gca,'show','Location','SouthEast')

% figure()
subplot(223)
hold on
for i = 1:length(r0)
    plot(sigma,meanRadius(i,:),'DisplayName',['\sigma = ',num2str(sigma(i))]);
end
hold off
xlabel('\sigma Disturbance Strength In'); ylabel('Mean Radius Out');
legend(gca,'show','Location','SouthEast')

% figure()
subplot(224)
hold on
for i = 1:length(r0)
    plot(sigma,stdRadius(i,:),'DisplayName',['r_0 = ',num2str(r0(i))])
end
hold off
xlabel('\sigma Disturbance Strength In'); ylabel('\sigma Radius Out');
legend(gca,'show','Location','SouthEast')

%% Plot the last realization

mRadius = meanRadius(end,end);
sRadius = stdRadius(end,end);

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
xc_m = mRadius*cos(tc)+xfoc_nom;
yc_m = mRadius*sin(tc)+yfoc_nom;
plot(xc_m,yc_m)
xc_std1 = (mRadius - sRadius)*cos(tc) + xfoc_nom;
yc_std1 = (mRadius - sRadius)*sin(tc) + yfoc_nom;
xc_std2 = (mRadius + sRadius)*cos(tc) + xfoc_nom;
yc_std2 = (mRadius + sRadius)*sin(tc) + yfoc_nom;
plot(xc_std1,yc_std1,'--k',xc_std2,yc_std2,'--k')
hold off
xlim([xfoc_nom-2*mRadius, xfoc_nom+2*mRadius]);
ylim([yfoc_nom-2*mRadius, yfoc_nom+2*mRadius]);
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
% Start = [x,y];
% Stop  = [xfoc,yfoc];
% 
% figure()
% arrow(Start,Stop,'Length',0.3,'BaseAngle',15,'TipAngle',25)
% axis square; axis equal

% %% Frequency analysis
% 
% dt = t(2);      % (s) sample time
% Fs = 1/dt;      % (Hz) sample frequecy
% Lx = length(x);
% NFFT = 2^nextpow2(Lx); % Next power of 2 from length of x
% X = fft(x,NFFT)/Lx;
% freq = Fs/2*linspace(0,1,NFFT/2+1);
% 
% Plot single-sided amplitude spectrum.
% figure()
% plot(freq,2*abs(X(1:NFFT/2+1))) 
% title('Single-Sided Amplitude Spectrum of y(t)')
% xlabel('Frequency (Hz)')
% ylabel('|Y(f)|')
