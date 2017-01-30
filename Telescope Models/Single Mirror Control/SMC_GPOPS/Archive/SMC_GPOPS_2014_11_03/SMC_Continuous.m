%---------------------------------------------%
% BEGIN: function SMC_Continuous.m %
%---------------------------------------------%
function phaseout = SMC_Continuous(input)
% the "output" from the continuous function includes "dynamics" (states)
% "path" (path constraints), and "integrand" (trajectory constraints).
% sprintf('Entering Continuous')

global integrand

%% Dynamics
% store auxiliary data locally
f = input.auxdata.f;
G = input.auxdata.G;
J = input.auxdata.J;
% w = auxdata.w;
ps = input.auxdata.ps;
% tw = auxdata.tw;
% continuous disturbance function
% fw = polyfit(w,tw,10);

% store time locally
t  = input.phase.time;

% store states locally
x1 = input.phase.state(:,1);
x2 = input.phase.state(:,2);
fu = input.phase.state(:,3);

% store controls locally
u = input.phase.control(:,1);

% dynamics equations
x1dot = x2;
% x2dot = (1/J)*(-G*x2 + fu + sin(100*t+ps) + sin(250*t-ps));     % let disturbance input be sin(omega*t)
x2dot = (1/J)*(-G*x2 + fu + 1*sin(10*t+ps));
fudot = u;
% fudot = zeros(size(u));

phaseout.dynamics = [x1dot, x2dot, fudot];

%% Path Constraint

% phaseout.path = vnorm;

%% Integral Constraint, Error Minimization
% There is no integral constraint
% integrand for minimizing spot size radius
% for i = 1:length(x1)
%     spotSize(i,1) = SingleMirror_SpotSize(x1(i),f);
% end
% integrand =  spotSize;
integrand = x1;
% size(integrand)
phaseout.integrand = 1e5*abs(integrand);

% sprintf('Leaving Continuous')
end
%---------------------------------------------%
% BEGIN: function SMC_Continuous.m %
%---------------------------------------------%
