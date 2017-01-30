%---------------------------------------------%
% BEGIN: function Q4_Continuous.m %
%---------------------------------------------%
function phaseout = Q2_Continuous(input)
% the "output" from the continuous function includes "dynamics" (states)
% "path" (path constraints), and "integrand" (trajectory constraints).
% sprintf('Entering Continuous')

%% Dynamics
% store time locally
t  = input.phase.time;
N = length(t);

% store states locally
x = input.phase.state(:,1);
y = input.phase.state(:,2);

% store controls locally
th = input.phase.control(:,1);      % sail angle
A =  input.phase.control(:,2);      % heading angle

% calculate absolute velocity vector of the boat
% velocities
% Vvector = zeros(2,N);
Vguess = [cos(A)'; sin(A)'];
Vriver = zeros(2,N);
 for i = 1:N
    if isnan(y(i))
        Vriver(1,i) = 0;
    else
        Vriver(1,i) = calcVriver(y(i));
    end
    if isnan(th(i)) || isnan(A(i))
        V = [0; 0];
        Vbr = [0; 0];
    else
        [V,Vbr] = calcVvector(th(i),A(i),Vguess(:,i),Vriver(:,i));
    end
    Vguess(:,i+1) = V;
    Vx(i,1) = V(1);
    Vy(i,1) = V(2);
end

% continuous "dynamics" equations
xdot = Vx;
ydot = Vy;

phaseout.dynamics = [xdot, ydot];

%% Path Constraint
% boat velocity wrt river must be in the direction of the heading angle
% (make sure the boat only goes forward)
% Vbr_angle = calcAngle(Vbr);
% Peq = [cos(A) - cos(Vbr_angle), sin(A) - sin(Vbr_angle)];
% phaseout.path = Peq;
% 
% % nested function to calculate the angle wrt +x-axis
% function angle = calcAngle(v)
%     % calc angle between v and +x-axis=
%     vx = v(1,:); vy = v(2,:);
%     angle = atan(vy./vx)';
%     angle(isnan(angle)) = 0;
% end

%% Integral Constraint, Error Minimization
% There is no integral constraint
% cost functional integrand

% sprintf('Leaving Continuous')
end
%---------------------------------------------%
% END: function Q4_Continuous.m %
%---------------------------------------------%
