function [ c, ceq ] = Constr_MECH622_Qual( u, N, s0 )
%CONSTR_MECH622_QUAL Contains all problem constraints
%   State Constraints
%   Initial Constraints
%   Endpoint Constraints
%   Path Constraints
%   c(x) <= 0, ceq(x) == 0

global w

% nsVars = min(size(s0));     % number of state variables

[s, Vbr] = States_MECH622_Qual(u, N, s0);

% % initial constraints
Ieq = s0(:,1) - s(:,1);

% endpoint constraints
Eeq = s0(:,end) - s(:,end);

% path constraints
% Add another constraint to prevent the force on the boat from
% being backwards. Velocity vector must be in the same direction as the
% heading angle at all times.
A = u(N+1:2*N);
Vbr_angle = calcAngle(Vbr);
Peq = [cos(A) - cos(Vbr_angle); sin(A) - sin(Vbr_angle)];

% bound the states to keep boat in the river
Py1 = -s(2,:);          % keep boat above lower bank
Py2 = s(2,:) - w;       % keep boat below upper bank

% assemble constraint vectors
ceq = [Ieq; Eeq; Peq];          % equality constraints vector

c = [Py1'; Py2'];               % inequality constraints vector

%% nested function to calculate the angle wrt +x-axis
function angle = calcAngle(v)
    % calc angle between v and +x-axis=
    vx = v(1,:); vy = v(2,:);
    angle = atan(vy./vx)';
    I = find(isnan(angle));
    angle(I) = 0;
end

end