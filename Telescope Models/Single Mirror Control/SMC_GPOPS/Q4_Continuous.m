%---------------------------------------------%
% BEGIN: function Q4_Continuous.m %
%---------------------------------------------%
function phaseout = Q4_Continuous(input)
% the "output" from the continuous function includes "dynamics" (states)
% "path" (path constraints), and "integrand" (trajectory constraints).
% sprintf('Entering Continuous')

% global integrand

%% Dynamics
% store x-position locally
x  = input.phase.time;

% store states locally
y = input.phase.state(:,1);

% store controls locally
u = input.phase.control(:,1);

% dynamics equations
ydot = u;

phaseout.dynamics = [ydot];

%% Path Constraint

% phaseout.path = vnorm;

%% Integral Constraint, Error Minimization
% There is no integral constraint
% integrand for energy (power)
% E = (1+u.^2).*(1+exp(-((x-2).^2+(y-2).^2)));
E = (1+exp(-((x-2).^2+(y-2).^2)));  % power w/o dir change cost

phaseout.integrand = E;

% sprintf('Leaving Continuous')
end
%---------------------------------------------%
% END: function Q4_Continuous.m %
%---------------------------------------------%
