% Tim Coon
% Qualifying Exam Question #2, MECH 622
% 01/26/2014
% Adapted from Bryson P4_3_5
% Solve the sailboat problem by an indirect numerical method wherein the
% control vector is guessed, then the states are propagated forward, the
% terminal constraints are used to find final costates, then the costates 
% are propagated backwards to find the derivative of the cost function wrt
% the control at each time step. These derivatives are used by the gradient
% search optimization algorithm of fmincon() to determine the search 
% direction that reduces the cost function. (Ref Bryson Sec 3.2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simplified Model - Formal Statement of the Optimization Problem
% Minimize:     J = tf
% State Eq:     x(i+1) = x(i) + Vx(i)*dt
%               y(i+1) = y(i) + Vy(i)*dt
%               Vx(i) = c1*Vwxr(i)*sin(th(i)) + Vriver(y(i))
%                   Vwxr(i) = Vwx - Vriver(y(i))
%                   Vwx = Vw*cos(Aw)
%                   Vwy = c2*Vwy*cos(th(i))
% Constraints:  x(0) = x0  y(0) = y0
%               x(N) = xf  y(N) = yf
%               0 <= y(i) <= w   i = 1,2,3,...,N
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%------------ given ------------------------
global c1 c2 xf yf Vmax w Vwx Vwy
Vwind = 6*0.514;            % (m/s)
Awind = deg2rad(130);       % (rad)
Vwx = Vwind*cos(Awind);     % (m/s)
Vwy = Vwind*sin(Awind);     % (m/s)
Vmax = 32*0.514;            % (m/s) max river velocity
% Vmax = Vr_ms(i);
w = 200;    % (m)
c1 = 0.7;   % coefficient for Vx
c2 = 0.7;   % coefficient for Vy
%---final states (posiion)
xf=115;     % (m)
yf=w;       % (m)
%---Set initial states [u x v y]
x0 = 85;    % (m)
y0 = 0;     % (m)
ns=2;  % number of state variables

%------------ guesses ------------------------
%---Set initial guess for final time
tf0=100;
%---Set number of time steps
N=100;
%---Set initial guess for control input (just try all zeros)
nc=1;               % number of control inputs
th0=(pi/2)*ones(N,1);     % theta is the sail angle wrt global +x
%---Set initial guess for states
s0=[linspace(x0,xf,N); linspace(y0,yf,N)];

%------------ nlp setup ------------------------
%---Set number of equality and inequality constraints
ncons=2;        % total number of constraints (neqcons+nineqcons)
neqcons=2;      % number of equality constraints
%---Define design vector, independent design variables in fmincon
% [theta deltaT]
u0=[th0; tf0/N];
% u0_data = load('u0_guess_N=100.mat');
% u0 = u0_data.u;

%---Define a vector holding all of the dimensions
dims=[N ns nc ncons];
name='dvdp';
% bounds on the input
lb = Awind-(3*pi/4); ub = Awind-(pi/2);
% Set options for fmincon
options = optimset('GradObj','on','Display','Iter','GradConstr','off','MaxIter',1e6,'MaxFunEvals',1e7);
u = fmincon(@(u0)Obj_MECH622_Qual(u0,s0,name,dims),u0,[],[],[],[],[],[],...
                        @(u0)Constr_MECH622_Qual(u0,s0,name,dims),options);

%% Retrieve state vector at optimum, organize for plotting
[f,g,s]=Obj_MECH622_Qual(u,s0,name,dims);
x=s(1,:);
y=s(2,:);
th=[u(1:N)]*180/pi;
tf=N*u(N+1);
t=tf*[0:1/N:1];

%% Plot
figure(1)
suptitle('Sailboat Min Time Solution (GradObj)')
subplot(4,6,[1 2 3 4, 7 8 9 10, 13 14 15 16, 19 20 21 22])
plot(x,y,'linewidth',2)
xlabel('x'); ylabel('y')
axis square
subplot(4,6,5:6)
plot(t,x,'linewidth',2)
xlabel('time'); ylabel('x');
subplot(4,6,11:12)
plot(t,y,'linewidth',2)
xlabel('time'); ylabel('y');
subplot(4,6,17:18)
plot(t(1:end-1),th,'linewidth',2)
xlabel('time'); ylabel('\theta');
ax = subplot(4,6,23:24);
set(ax,'visible','off')
t_f = strcat('t_f = ', num2str(tf,'%6.4f'));
ht = text(0.2,0.3,t_f);
ht.FontSize = 16;

