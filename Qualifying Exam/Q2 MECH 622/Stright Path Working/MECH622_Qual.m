% Tim Coon
% Qualifying Exam Question #2, MECH 622
% 11/20/2014
% straightforward method

% clear all; 
close all; clc;

%------------ given ------------------------
global Vmax w Vw
Vwind = 6*0.514;            % (m/s) wind speed
Awind = deg2rad(130);       % (rad) wind direction angle wrt x-axis
Vw = [Vwind*cos(Awind); Vwind*sin(Awind)];
Vmax = 32*0.514;            % (m/s) max river velocity
w = 200;                    % (m) river width
%---final states (posiion)
xf=115;     % (m)
yf=w;       % (m)
%---initial states
x1 = 85;    % (m)
y1 = 0;     % (m)

%------------ guesses ------------------------
%---Set initial guess for final time
tf0=150;
%---Set number of time states
N=10;
%---Set initial guess for control input
nc=2;                % number of control inputs
th0=(-pi/3)*ones(N,1);     % theta is the sail angle wrt -boat centerline
A0 = (atan(20/3))*ones(N,1);     % A is the heading angle wrt x-axis
% A0 = linspace(0,pi/2,N)';
u0 = [th0; A0; tf0/N];

% initial guess for states
x0 = linspace(x1,xf,N); y0 = linspace(y1,yf,N);
s0 = [x0; y0];

%------------ nlp setup ------------------------

optn=optimset('GradConst','off','Display','Iter','MaxIter',100);
lb = zeros(length(u0),1); ub = zeros(length(u0),1);
lb(1:N) = (-pi/2)*ones(N,1); ub(1:N) = (pi/2)*ones(N,1);    % sail angle
lb(N+1:2*N) = zeros(1,N); ub(N+1:2*N) = 2*pi;               % heading angle
lb(2*N+1) = 100/N; ub(2*N+1) = 600/N;                       % time step
[u fval] = fmincon(@Cost_MECH622_Qual,u0,[],[],[],[],lb,ub,...
                                       @Constr_MECH622_Qual,optn,N,s0);
% u = u0;
th = u(1:N);
A = u(N+1:2*N);
dt = u(2*N+1);
t = 0:dt:(N-1)*dt;
tc = 0:dt:(N-1)*dt;

s = States_MECH622_Qual(u, N, s0);
x = s(1,:);
y = s(2,:);

%% Plot results
figure(1)
plot(x,y,'linewidth',2)
xlabel('x'); ylabel('y')

figure(2)
subplot(211)
plot(t,x,'linewidth',2)
xlabel('time'); ylabel('x');
subplot(212)
plot(t,y,'linewidth',2)
xlabel('time'); ylabel('y');

figure(3)
subplot(211)
plot(tc,th,'linewidth',2)
xlabel('time'); ylabel('\theta');
subplot(212)
plot(tc,A,'linewidth',2)
xlabel('time'); ylabel('A');
