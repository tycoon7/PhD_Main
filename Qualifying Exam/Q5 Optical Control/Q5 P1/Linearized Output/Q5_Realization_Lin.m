function [t,z,Dt,D,X] = Q5_Realization_Lin(r0,sigma0)
% generate a random disturbance and dynamic response realization as well as
% the linerized output
% r0     = disturbance strength (mean)
% sigma0 = disturbance standard deviation
% C_WF   = linearized output matrix per RedBreck
% t      = dynamics time vector
% z      = first-order state vector
% Dt     = disturbance time vector
% D      = disturbance vector
% Xn     = output vector (all of the output rays in form of Eq 113) 

global t0 tf

% Initial Conditions
X0  = 0;    xd0 = 0;
y0  = 0;    yd0 = 0;
a0  = 0;    ad0 = 0;
z0 = [X0 xd0 y0 yd0 a0 ad0]';

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
tspan   = linspace(t0,tf,100);    % this vector controls the number of states
options = odeset('RelTol',1e-2);  % set tolerances for integration
[t,z]   = ode45(@state_eqns,tspan,z0,options,Dt,D);
x       = z(:,1);    % (m) CoM x-position
% xd     = z(:,2);    % (m/s) CoM x-velocity
y       = z(:,3);    % (m) CoM y-position
% yd     = z(:,4);    % (m/s) CoM y-velocity
alpha   = z(:,5);    % (rad) mirror angular displacement
% alphad = z(:,6);    % (rad/s) mirror angular rate

%% Linearized Output Equations
% Arrange the states to work with the wavefront model of RedBreck Eq. 121.
% Each column is one state.
% the input beam does not change in this simulation
n = length(t);      % n states
X0 = zeros(7,n);
% the dynamics of surface #1 are arranged from the states above
theta1 = [zeros(1,n); zeros(1,n); alpha'];
delta1 = [x'; y'; zeros(1,n)];
u1 = [theta1; delta1];
% surface #2 does not change in this simulation
u2 = zeros(6,n);
% the complete state vectors for output calcs are
X = [X0; u1; u2];
% the output dynamics are easily calculated
% Xn = C_WF*X;
