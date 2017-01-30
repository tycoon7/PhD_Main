% Tim Coon - 27 February, 2015
% Qualifying Exam, Question #5, Part 1
% With no active control, characterize the evolution of the focal point
% distribution in response to the random disturbance input.
clear; close all; clc;

global f_par rho1 rho2 phi1 phi2 Pc phi Pc_mag P0 I0 I m c k t0 tf

% Geometry
x1 = -0.4;                  % (m) left edge of mirror
D = 0.25;                   % (m) clear aperture of mirror
x2 = x1 + D;                % (m) right edge of mirror
f_par = 1;                  % (m) focal length of mirror
den = 0.1;                  % (kg/m) length density of mirror
[Pc,P1,P2,m,I] = calcConicCentroid(x1,x2,f_par,den);
rho1 = P1 - Pc;                     % (m) vector from CoM to point 1
rho2 = P2 - Pc;                     % (m) vector from CoM to point 2
phi1 = atan(abs(rho1(1)/rho1(2)));     % (m) angle of rho1 wrt  y-axis
phi2 = atan(abs(rho2(1)/rho2(2)));     % (m) angle of rho2 wrt -y-axis
phi  = atan(abs(Pc(1)/Pc(2)));         % {m) angle of Pc wrt -y-axis
Pc_mag = norm(Pc);                 % (m) magnitude of CoM vector wrt vertex
cos_phi1 = cos(phi1);
sin_phi1 = sin(phi1);
cos_phi2 = cos(phi2);
sin_phi2 = sin(phi2);
P0 = [x1+D/2; 4; 0];       % starting point of chief ray
I0 = [0; -1; 0];            % initial direction of chief ray

% Mirror Parameters
f1 = f_par;                         % (m) focal length of parabolic mirror
f2 = inf;                           % (m) focal length of reference surface
e1 = 1;                             % (-) eccentricity of par
e2 = 0;                             % (-) eccentricity of rs
V1 = [0.00; 0.00;  0.00];         % (m) vertex, par
V2 = [0.00;  f1;  0.00];            % (m) vertex, rs
q1 = [Pc(1); Pc(2); 0.00];       % (m) rotation point, par
q2 = V2;                            % (m) rotation point, rs
psi1 = [0;  1; 0];                  % (-) principal axis of par
psi2 = [0; -1; 0];                  % (-) principal axis of rs

% Assemble parameters
e = [e1 e2];
f = [f1 f2];
V = [V1 V2];
q = [q1 q2];
psi = [psi1 psi2];

P0_chief = [(x1+x2)/2; 1.5; 0];
I0 = [0; -1; 0];

% Dynamics Properties
% How do I "scale" these?
c = 0.1;
k = 5;
t0 = 0;
tf = 1;

%% build mirror objects
for s = 1:length(e)
    M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),eye(3));
%     M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),rotz(-Thz(s)));
end

%% build beam objects (just one)
rd = 3;         % ray density
beam = BeamClass(D,rd,P0_chief,I0,M);

%% generate sample realizations and determine statistics of resulting focal-
% point-evolution ball radius
nreals = 7;
radius = zeros(nreals,1);
r0 = 0.1:0.1:0.3;
sigma0 = 0.05:0.05:0.15;
meanAvgSpotSize  = zeros(length(r0),length(sigma0));
stdAvgSpotSize  = zeros(length(r0),length(sigma0));

%% simulate nreals realizations for each combination of r0 and sigma0
for i = 1:length(r0)
    for j = 1:length(sigma0)
        for n = 1:nreals
            tic
            [t,z,Dt,D,X] = Q5_Realization_Lin(r0(i),sigma0(j));
            linearWFmodel(X,beam);
            toc
        end
        meanAvgSpotSize(i,j) = mean(mean(spotSize,1));
        stdAvgSpotSize(i,j) = std(mean(spotSize,1));
    end
end

x      = z(:,1);
xd     = z(:,2);
y      = z(:,3);
yd     = z(:,4);
alpha  = z(:,5);
alphad = z(:,6);
Dx     = D(:,1);
Dy     = D(:,2);

save('sim_7reals_100states_randImpulse.mat')


%% Plot out stats wrt in stats
set(0,'DefaultAxesFontSize', 20)            % make default plot font bigger
% set(0,'defaultaxescolororder',[0 0 0; 0.5 0.5 0.5]) %black and gray
% set(0,'defaultaxeslinestyleorder',{'-','--','-*',':'}) %or whatever you want
load('sim_7reals_100states_randImpulse.mat')
mrksData = ['- ';'--';'-*';': '];
mrks = cellstr(mrksData);
% mu vs r0
statPlot(1) = figure();
% suptitle('Output Stats Vs. Input Stats')
hold on
for i = 1:length(sigma0)
    plot(r0,meanAvgSpotSize(:,i),char(mrks(i)),'DisplayName',['\sigma_0 = ',num2str(sigma0(i))],'Linewidth',4);
end
hold off
xlabel('Mean Disturbance Strength In'); ylabel('Mean Spot Size Out');
legend(gca,'show','Location','EastOutside')

% sigma vs r0
statPlot(2) = figure();
hold on
for i = 1:length(sigma0)
    plot(r0,stdAvgSpotSize(:,i),char(mrks(i)),'DisplayName',['\sigma_0 = ',num2str(sigma0(i))],'Linewidth',4)
end
hold off
xlabel('Mean Disturbance Strength In'); ylabel('Std Spot Size Out');
legend(gca,'show','Location','EastOutside')

% mu vs sigma0
statPlot(3) = figure();
hold on
for i = 1:length(r0)
    plot(sigma0,meanAvgSpotSize(i,:),char(mrks(i)),'DisplayName',['r_0 = ',num2str(r0(i))],'Linewidth',4);
end
hold off
xlabel('Std Disturbance Strength In'); ylabel('Mean Spot Size Out');
legend(gca,'show','Location','EastOutside')

% sigma vs sigma0
statPlot(4) = figure();
hold on
for i = 1:length(r0)
    plot(sigma0,stdAvgSpotSize(i,:),char(mrks(i)),'DisplayName',['r_0 = ',num2str(r0(i))],'Linewidth',4)
end
hold off
xlabel('Std Disturbance Strength In'); ylabel('Std Spot Size Out');
legend(gca,'show','Location','EastOutside')

savefig(statPlot,'StatFigs_7reals_100states_randImpulse.fig')     % open with openfig()
print(statPlot(1),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_stat_mu_v_r0.pdf')
print(statPlot(2),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_stat_sigma_v_r0.pdf')
print(statPlot(3),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_stat_mu_v_sigma0.pdf')
print(statPlot(4),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_stat_sigma_v_sigma0.pdf')

%% Plot the last realization

% plot CoM motion
realPlot(1) = figure();
plot(x,y)
xlabel('x-pos (m)'); ylabel('y-pos (m)');
axis equal

% plot CoM directional displacements
realPlot(2) = figure();
plot(t,x,'b-',t,y,'r--','Linewidth',4)
legend('x','y')
xlabel('time'); ylabel('Force');

% plot mirror angular displacement
realPlot(3) = figure();
plot(t,alpha,'Linewidth',4)
xlabel('time (s)'); ylabel('Angle (rad)');

% plot the disturbance input
realPlot(4) = figure();
plot(Dt,Dx,'b-',Dt,Dy,'r--','Linewidth',4)
legend('Dx','Dy')
xlim([-0.1, tf]);
xlabel('time'); ylabel('Displacement');

% plot the evolution of the spot size and the average
realPlot(5) = figure();
plot(t,spotSize(:,1),'b-',t,mean(spotSize(:,1),1)*ones(length(t),1),'r--','Linewidth',4)
legend('Spot Size','Avg Size')
xlabel('time'); ylabel('Spot Size');

savefig(realPlot,'RealFigs_7reals_100states_randImpulse.fig')     % open with openfig()

print(realPlot(1),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_real_CoM_Motion.pdf')
print(realPlot(2),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_real_CoM_LinDisp.pdf')
print(realPlot(3),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_real_CoM_AngDisp.pdf')
print(realPlot(4),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_real_Disturbance.pdf')
print(realPlot(5),'-dpdf','.\Q5P1_LaTeX\Figures\Q5P1_real_SpotSizeEvolution.pdf')