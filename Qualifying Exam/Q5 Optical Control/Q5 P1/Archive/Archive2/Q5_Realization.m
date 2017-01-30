function [t,z,radius,spotSize,D] = Q5_Realization()
% generate a random disturbance and dynamics realization

global f Pc phi Pc_mag r0 sigma t0 tf

% Initial Conditions
x0  = 0;    xd0 = 0;
y0  = 0;    yd0 = 0;
a0  = 0;    ad0 = 0;
z0 = [x0 xd0 y0 yd0 a0 ad0]';

%% Create a new disturbance input for each realization
Dt = linspace(t0,tf)';
% generate random variables
r = sigma*randn(length(Dt),1) + r0;
theta = 2*pi*rand(length(Dt),1);
% random impulse input vector
D = zeros(length(Dt),2);
D(1,:) = [r(1)*cos(theta(1)), r(1)*sin(theta(1))];
% random continuous input vector
% D = [r.*cos(theta), r.*sin(theta)];
% continuous sinusoidal input vector
% D = [sin(1*Dt), sin(2*Dt)];

%% integrate state equations
tspan = linspace(t0,tf);
options=odeset('RelTol',1e-2);  % set tolerances for integration
[t,z] = ode45(@state_eqns,tspan,z0,options,Dt,D);
x      = z(:,1);    % CoM x-position
xd     = z(:,2);    % CoM x-velocity
y      = z(:,3);    % CoM y-position
yd     = z(:,4);    % CoM y-velocity
alpha  = z(:,5);    % mirror angular displacement
alphad = z(:,6);    % mirror angular rate

%% calculate the focal point positions by first finding the vertex position
% from the CoM, then find focal point from vertex 
% mirror vertex position
Vx_nom = -Pc(1);
Vy_nom = -Pc(2);
Vx = Vx_nom + x + Pc_mag*alpha*cos(phi);
Vy = Vy_nom + y + Pc_mag*alpha*sin(phi);
dV = [Vx-Vx_nom Vy-Vy_nom];

% focal point position
xfoc_nom = Vx_nom;
yfoc_nom = Vy_nom + f;
xfoc = Vx + f*sin(alpha);
yfoc = Vy + f*cos(alpha);

%% Calculate radius of focal point evolution ball centered on the nominal
d = zeros(length(t),1);
for i = 1:length(t)
    d(i) = norm([xfoc_nom - xfoc(i), yfoc_nom - yfoc(i)]);
end

radius = max(d);

%% Calculate the spot size at each state
for c = 1:length(t)
    spotSize(c,1) = SingleMirror_SpotSize(alpha(c),dV(c,:));
end

disp(spotSize)