function [t,z,spotSize,D] = zQ5_P1_Realization(r0,sigma0)
% generate a random disturbance and dynamics realization

global Pc phi Pc_mag t0 tf

% Initial Conditions
x0  = 0;    xd0 = 0;
y0  = 0;    yd0 = 0;
a0  = 0;    ad0 = 0;
z0 = [x0 xd0 y0 yd0 a0 ad0]';

%% Create a new disturbance input for each realization
Dt = linspace(t0,tf)';
% generate random variables
r = sigma0*randn(length(Dt),1) + r0;
theta = 2*pi*rand(length(Dt),1);
% random impulse input vector
D = zeros(length(Dt),2);
D(1,:) = [r(1)*cos(theta(1)), r(1)*sin(theta(1))];
% random continuous input vector
% D = [r.*cos(theta), r.*sin(theta)];
% continuous sinusoidal input vector
% D = [sin(1*Dt), sin(2*Dt)];

%% integrate state equations
tspan = linspace(t0,tf,100);    % this vector controls the number of states
options=odeset('RelTol',1e-2);  % set tolerances for integration
[t,z] = ode45(@state_eqns,tspan,z0,options,Dt,D);
x      = z(:,1);    % (m) CoM x-position
% xd     = z(:,2);    % (m/s) CoM x-velocity
y      = z(:,3);    % (m) CoM y-position
% yd     = z(:,4);    % (m/s) CoM y-velocity
alpha  = z(:,5);    % (rad) mirror angular displacement
% alphad = z(:,6);    % (rad/s) mirror angular rate

%% mirror vertex position
Vx_nom = -Pc(1);
Vy_nom = -Pc(2);
Vx = Vx_nom + x + Pc_mag*alpha*cos(phi);
Vy = Vy_nom + y + Pc_mag*alpha*sin(phi);
dV = [Vx-Vx_nom Vy-Vy_nom];

%% Calculate the Strehl Ratio at each state using Zemax
