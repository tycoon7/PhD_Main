function [s, Vbr_hat] = States_MECH622_Qual(u, N, s0)
% STATES_MECH622_QUAL calculates the states

global Vw mb
x = s0(1,:);
y = s0(2,:);
V = s0(3:4,:);
th = u(1:N);
A = u(N+1:2*N);
dt = u(2*N+1);
c1 = 0.5;     % wind force coeff
c2 = 0.5;     % water drag force coeff

% velocities
Vdot = zeros(2,N);
Vbr_hat = zeros(2,N);
Vriver = zeros(2,N);

% state propagation
for i = 1:N-1
    Vriver(1,i) = calcVriver(y(i));
    [Vdot(:,i), Vbr_hat(:,i)] = calcVdot(A(i),th(i),V(:,i),Vriver(:,i));
    x(i+1) = x(i) + V(1,i)*dt;
    y(i+1) = y(i) + V(2,i)*dt;
    V(:,i+1) = V(:,i) + Vdot(:,i)*dt;
end

s = [x; y; V];

%% Nested function to calculate force on the boat
function [ Vdi, Vbr_hati] = calcVdot(Ai,thi,Vi,Vri)
    Vwbi = Vw - Vi;
    Vbri = Vi - Vri;
    rs = [cos((Ai+pi)+thi); sin((Ai+pi)+thi)];    % main sail direction vector
    Vbr_hati = [cos(Ai); sin(Ai)];
    rs_X_Vwb = det([rs';Vwbi']);
    sine_theta = det([-Vbr_hati';rs']);
    Fs_mag = c1*norm(Vwbi)*rs_X_Vwb;
    Fs_br = Fs_mag*sine_theta*Vbr_hati;
    Fd = c2*norm(Vbri)*Vbri;
    Fb = Fs_br - Fd;
    Vdi = (1/mb)*Fb;
end

end