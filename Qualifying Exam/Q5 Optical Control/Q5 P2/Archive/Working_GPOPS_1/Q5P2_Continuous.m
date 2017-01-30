%---------------------------------------------%
% BEGIN: function Q5P2_Continuous.m %
%---------------------------------------------%
function phaseout = Q5P2_Continuous(input)
% the "output" from the continuous function includes "dynamics" (states)
% "path" (path constraints), and "integrand" (trajectory constraints).
% sprintf('Entering Continuous')

%% Dynamics
% store auxiliary data locally
% R = input.auxdata.R;
R1 = input.auxdata.R1;
R2 = input.auxdata.R2;
rho = input.auxdata.rho;
I = input.auxdata.I;
m = input.auxdata.m;
c = input.auxdata.c;
k = input.auxdata.k;
r0 = input.auxdata.r0;
sigma = input.auxdata.sigma;
cos_phi1 = input.auxdata.cos_phi1;
sin_phi1 = input.auxdata.sin_phi1;
cos_phi2 = input.auxdata.cos_phi2;
sin_phi2 = input.auxdata.sin_phi2;

% store time locally
t  = input.phase.time;

% store states locally
x = input.phase.state(:,1);
y = input.phase.state(:,2);
alpha = input.phase.state(:,3);
xd = input.phase.state(:,4);
yd = input.phase.state(:,5);
alphad = input.phase.state(:,6);

% store controls locally
Ux = input.phase.control(:,1);
Uy = input.phase.control(:,2);

% disturbance
Dx = sin(t);
Dy = sin(t);

% dynamics equations

x1  = x  - rho*alpha*cos_phi1;
x1d = xd - rho*alphad*cos_phi1;
x2  = x  + rho*alpha*cos_phi2;
x2d = xd + rho*alphad*cos_phi2;
y1  = y  - rho*alpha*sin_phi1;
y1d = yd - rho*alphad*sin_phi1;
y2  = y  + rho*alpha*sin_phi2;
y2d = yd + rho*alphad*sin_phi2;
x0  = x1 + Ux;
x0d = x1d;
y0  = y1 + Uy;
y0d = y1d;

xdot = xd;
xdd = (1/m)*(-c*x0d - k*x0 - c*x2d - k*x2 + Dx);
ydot = yd;
ydd = (1/m)*(-c*y0d - k*y0 - c*y2d - k*y2 + Dy);
adot = alphad;
add = (1/I)*(R1*(c*x0d + k*x0 - c*y2d - k*y2 + Dy)...
               + R2*(-c*x2d - k*x2 + Dx + c*y0d + k*y0));

phaseout.dynamics = [xdot, xdd, ydot, ydd, adot, add];

%% Path Constraint

% phaseout.path = vnorm;

%% Integral Constraint, Error Minimization
% There is no integral constraint
% integrand for minimizing spot size radius
% for i = 1:length(x1)
%     spotSize(i,1) = SingleMirror_SpotSize(x1(i),f);
% end
% integrand =  spotSize;
integrand = sqrt(x.^2+y.^2)+alpha.^2;
% size(integrand)
phaseout.integrand = abs(integrand);

% sprintf('Leaving Continuous')
end
%---------------------------------------------%
% END: function Q5P2_Continuous.m %
%---------------------------------------------%
