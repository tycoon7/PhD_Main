% Tim Coon
% Qualifying Exam Question #2, MECH 622
% 01/12/2014
% straightforward method

clear; close all; clc;

%------------ given ------------------------
global Vmax w Vw mb
Vwind = 6*0.514;            % (m/s) wind speed
Awind = deg2rad(130);       % (rad) wind direction angle wrt x-axis
Vw = [Vwind*cos(Awind); Vwind*sin(Awind)];
Vmax = 32*0.514;            % (m/s) max river velocity
w = 200;                    % (m) river width
mb = 200;                  % (kg) mass of boat
%---final states (posiion)
xf=115;         % (m)
yf=w;           % (m)
%---initial states
xi = 85;        % (m)
yi = 0;         % (m)
Vi = [0; 0];    % (m/s)

%------------ guesses ------------------------
%---Set initial guess for final time
tf0=200;
%---Set number of time steps
N=100;
%---Set initial guess for control input
nc=2;                           % number of control inputs
th0=(pi/2)*ones(N,1);          % theta is the sail angle wrt -boat C/L
A0 = (atan(20/3))*ones(N,1);   % A is the heading angle wrt x-axis
u0 = [th0; A0; tf0/N];

% initial guess for states
x0 = linspace(xi,xf,N); y0 = linspace(yi,yf,N);
V0 = ones(2,N); V0(:,1) = Vi;
s0 = [x0; y0; V0];

%------------ nlp setup ------------------------

optn=optimset('GradConst','off','Display','Iter','MaxIter',1e6,'MaxFunEvals',1e6,'algorithm','sqp');
lb = zeros(length(u0),1); ub = zeros(length(u0),1);
lb(1:N) = -pi/2; ub(1:N) = pi/2;                % sail angle
lb(N+1:2*N) = -2*pi; ub(N+1:2*N) = 2*pi;        % heading angle
lb(2*N+1) = 10/N; ub(2*N+1) = 1000/N;            % time step
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
plot(tc,th,'linewidth',2)
xlabel('time'); ylabel('\theta');
subplot(4,6,23:24)
plot(tc,A,'linewidth',2)
xlabel('time'); ylabel('A');