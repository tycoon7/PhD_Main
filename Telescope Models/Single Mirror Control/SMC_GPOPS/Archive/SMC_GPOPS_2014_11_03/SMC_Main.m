%---------------------------------------------------%
% Single Mirror Control:                            %
% AFIT AERO 899 Summer 2014 Dr Cobb                 %
% Timothy Coon                                      %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
% (see Burl Ch 2, C.E. 2.1)                         %
%   Minimize spot size                              %
% subject to the dynamic constraints                %
%    dx1/dt = x2                                    %
%    dx2/dt = (1/J)*(-G*x2 + u + w)                 %   
%    x1 = angle error                               %
%    x2 = angle error rate                          %
% and the boundary conditions                       %
%    x1(0) = 0, x2(0) = 0                           %
%    x1(t_f) = FREE, x2(t_f) = FREE                 %
%---------------------------------------------------%

clear all;
close all; clc;

auxdata.f = 0.70;              % (m) focal length of parabolic mirror
auxdata.G = 0.05;             % (kg-m) friction force per (m/s)
auxdata.J = 4*0.2/(3*pi);   % (kg-m^2) MOI of primary mirror
% auxdata.w = 3.*randn(100,1);     % (N) random disturbance force
% auxdata.ps = pi*rand;         % (rad) random phase shift
auxdata.ps = 0;
                                
% time boundary conditions
t0 = 0;                       % (s) force start at time = zero
tf = 1;                       % (s) fixed simulation time
% auxdata.tw = linspace(t0,tf,length(auxdata.w));

% state limits
x1min = -2; x1max = 2;        % (m) min/max angle error
x2min = -5; x2max = 5;        % (m) min/max angle error rate
fumin = -10; fumax = 10;        % (N) min/max control force

% control limits
% umin = -1; umax = 1;          % (N-m) control torque
umin = -10; umax = 10;          % (W) control power

% state boundary conditions
x10 = 0;                      % (m)
x20 = 0;                      % (m)
fu0 = 0;                      % (N)

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
iphase = 1;
% time
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;

% states
% fixed initial state
bounds.phase.initialstate.lower = [x10,x20,fu0]; 
bounds.phase.initialstate.upper = [x10,x20,fu0];
% bounds for states after initial
bounds.phase.state.lower = [x1min,x2min,fumin]; 
bounds.phase.state.upper = [x2max,x2max,fumax];
% free final state
bounds.phase.finalstate.lower = [x1min,x2min,fumin]; 
bounds.phase.finalstate.upper = [x1max,x2max,fumax];

% control limits
bounds.phase.control.lower = [umin]; 
bounds.phase.control.upper = [umax];

% path constraints
% bounds.phase.path.lower = 1;
% bounds.phase.path.upper = 1;
% The integral constraints are inequalities. Express this by setting the
% lower bound to be the minimum required, then a guess for the upper
% bound which is NOT ARBITRARILY HIGH, MUST BE "CLOSE". 
bounds.phase.integral.lower = [0];
bounds.phase.integral.upper = [100];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t0; tf];
guess.phase.state   = [[x10, x10, fu0];[x20, x20, fu0]];
guess.phase.control = [[0; 0]];
% the integral guess is to be a row vector (1xnd)
guess.phase.integral = [0];

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name = 'SMC-Problem';
setup.functions.continuous = @SMC_Continuous;
setup.functions.endpoint = @SMC_Endpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;
setup.nlp.solver = 'ipopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'first';
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-3;
% setup.mesh.maxiteration = 2;
% setup.mesh.colpointsmin = 4;
% setup.mesh.colpointsmax = 4;
setup.method = 'RPMintegration';
% setup.mesh.phase.colpoints = 4*ones(1,10);
% setup.mesh.phase.fraction = 0.1*ones(1,10);
% setup.scales.method = 'automatic-hybridUpdate';

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
output.result.nlptime
solution = output.result.solution;
GPOPS2_ExitFlags(output)

%% -----------------------------------------------------------------------%
%------------------------------- Plot Solution ---------------------------%
%-------------------------------------------------------------------------%
figure('position',[100 50 1100 800])
subplot(221)
% plot the states
p1 = plot(solution.phase(1).time, solution.phase(1).state(:,1:2), '-o');
xl = xlabel('time');
yl = ylabel('state');
set(p1,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
h1 = legend('$e$','$\dot{e}$','Location','eastoutside');   
set(h1,'interpreter','latex')
grid on

subplot(222)
% plot the control
p2 = plot(solution.phase(1).time,solution.phase(1).control,'-o');
xl = xlabel('time');
yl = ylabel('control');
set(p2,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
h2 = legend('$u$','Location','eastoutside');
set(h2,'interpreter','latex')
grid on

ps = auxdata.ps;
% w = sin(100*linspace(t0,tf,1000)+ps) + sin(250*linspace(t0,tf,1000)-ps);
w = 1*sin(10*linspace(t0,tf,1000)+ps);

subplot(223)
% plot the disturbance
% pp = plot(solution.phase(1).time,sin(10*solution.phase(1).time),'-o');
p3 = plot(linspace(t0,tf,1000),w,solution.phase.time,solution.phase.state(:,3));
xl = xlabel('time');
yl = ylabel('Force');
set(p3,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
h3 = legend('$w$','$F_u$','Location','eastoutside');
set(h3,'interpreter','latex')
grid on

% calculate the spot size wrt time using the final states
e = solution.phase.state(:,1);
spotSize = zeros(length(e),1);
for i = 1:length(e)
    spotSize(i,1) = SingleMirror_SpotSize(e(i),auxdata.f);
end
max_ss = max(spotSize);

hp4 = subplot(224);
focalLength = num2str(auxdata.f);
fString = strcat('Focal Length = ',focalLength);
text(0.25,0.8,fString);
inttext='Integral = ';
int = num2str(solution.phase.integral);
intstring = strcat(inttext,int);
text(0.25,0.6,intstring);
maxSStext = 'Maximum Spot Size = ';
maxSS = num2str(max_ss);
SSstring = strcat(maxSStext,maxSS);
text(0.25,0.4,SSstring);
set(hp4,'visible','off');

% % plot the cost integrand wrt time
% p4 = plot(solution.phase.time,spotSize);
% xl = xlabel('time');
% yl = ylabel('RMS Spot Size');
% % yl = ylabel('Angle Error');
% set(p4,'LineWidth',1.25,'MarkerSize',8);
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% % h4 = legend('$\tau$','Location','eastoutside');
% % set(h4,'interpreter','latex')
% grid on

filename = strcat('f_',focalLength,'.fig');
savefig(filename);
