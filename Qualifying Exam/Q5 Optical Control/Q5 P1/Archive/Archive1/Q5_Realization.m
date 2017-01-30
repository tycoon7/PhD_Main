function [t,z,radius,D] = Q5_Realization()
% generate a random disturbance and dynamics realization

global f CoM_dist rho2 phi2 r0 sigma t0 tf 

% Initial Conditions
x0  = 0;    xd0 = 0;
y0  = 0;    yd0 = 0;
a0  = 0;    ad0 = 0;
z0 = [x0 xd0 y0 yd0 a0 ad0]';

%% Create a new disturbance input for each realization
Dt = linspace(t0,tf,100)';
% generate random variables
r = sigma*randn(length(Dt),1) + r0;
theta = 2*pi*rand(length(Dt),1);
% random impulse input vector
D = zeros(length(Dt),2);
D(1,:) = [r(1)*cos(theta(1)), r(1)*sin(theta(1))];
% random continuous input vector
% D = [r.*cos(theta), r.*sin(theta)];
% continuous sinusoidal input vector
% D = [sin(10*Dt), sin(5*Dt)];

%% integrate state equations
options=odeset('RelTol',1e-3);  % set tolerances for integration
[t,z] = ode45(@state_eqns,[t0 tf],z0,options,Dt,D);
x      = z(:,1);    % CoM x-position
xd     = z(:,2);    % CoM x-velocity
y      = z(:,3);    % CoM y-position
yd     = z(:,4);    % CoM y-velocity
alpha  = z(:,5);    % mirror angular displacement
alphad = z(:,6);    % mirror angular rate

%% calculate the focal point positions by first finding the vertex position
% from the CoM, then find focal point from vertex 
% mirror vertex position
x2_nom = rho2(1);
y2_nom = rho2(2);
x2  = x2_nom + x  + CoM_dist*alpha*cos(phi2);
y2  = y2_nom + y  + CoM_dist*alpha*sin(phi2);
% focal point position
xfoc_nom = x2_nom;
yfoc_nom = y2_nom + f;
xfoc = x2 + f*sin(alpha);
yfoc = y2 + f*cos(alpha);

%% Calculate radius of focal point evolution ball centered on the nominal
d = zeros(length(t),1);
for i = 1:length(t)
    d(i) = norm([xfoc_nom - xfoc(i), yfoc_nom - yfoc(i)]);
end

radius = max(d);

end