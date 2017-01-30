% Tim Coon
% Qualifying Exam Question #2, MECH 622
% 11/17/2014
% Adapted from Bryson P4_3_5

clear all; close all; clc;

%------------ given ------------------------
global xf yf Vmax w Vwx Vwy
Vwind = 6*0.514;            % (m/s)
Awind = deg2rad(130);       % (rad)
Vwx = Vwind*cos(Awind);     % (m/s)
Vwy = Vwind*sin(Awind);     % (m/s)
Vmax = 32*0.514;            % (m/s) max river velocity
w = 200;    % (m)
%---final states (posiion)
xf=115;     % (m)
yf=w;       % (m)
%---Set initial states [u x v y]
x0 = 85;    % (m)
y0 = 0;     % (m)
s0=[x0 y0]';
ns=length(s0);

%------------ guesses ------------------------
%---Set initial guess for final time
tf=100;
%---Set number of time steps
N=200;
%---Set initial guess for control input (just try all zeros)
nc=2;               % number of control inputs
th=zeros(N,1)';     % theta is the sail angle wrt boat centerline
A = ones(N,1)';     % A is the heading angle wrt x-axis 

%------------ nlp setup ------------------------
%---Set number of equality and inequality constraints
ncons=2;        % total number of constraints (?)
neqcons=2;      % number of equality constraints
%---Define design vector, independent design variables in fmincon
% [theta deltaT]
thAt0=[th A tf/N];
%---Define a vector holding all of the dimenstions
dims=[N ns nc ncons];
% Name of evaluation function used by fcn and grd routines, we'll stick
% with Bryson's here
name='dvdp';
% Set options for fmincon
% options=optimset('Display','iter','GradObj','on','GradConstr','on','TolFun',1e-9,...
%     'TolCon',1e-9,'TolX',1e-5);
% With finite difference gradients
options = optimset('GradObj','off','GradConstr','off');
thAt=fmincon(@(thAt0)Obj_MECH622_Qual(thAt0,s0,name,dims),thAt0,[],[],[],[],[],[],@(thAt0)Constr_MECH622_Qual(thAt0,s0,name,dims),options);
%---Retrieve state vector at optimum, organize for plotting
[f,g,s]=Obj_MECH622_Qual(thAt,s0,name,dims);
u=s(1,:);
x=s(2,:);
v=s(3,:);
y=s(4,:);
thh=[thAt(1:N) thAt(N)]*180/pi;
tf=N*thAt(1,N+1);
t=tf*[0:1/N:1];

% % Plot
% figure(1);
% clf;
% plot(x,y,x,y,'b.',x(N+1),y(N+1),'ro');
% grid;
% % figure(1); clf; plot(x,y,x(N+1),y(N+1),'ro'); grid;
% xlabel('x');
% ylabel('y');
% %
% figure(2);
% clf;
% subplot(211), Bryson_ZOH_plot(t,thh);
% grid;
% axis([0 1 -50 0]);
% ylabel('\theta (deg)');
% subplot(212), plot(t,u,t,v);
% grid;
% xlabel('t');
% ylabel('u and v');

