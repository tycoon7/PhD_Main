function [s, Vvector] = States_MECH622_Qual(u, N, s0)
% STATES_MECH622_QUAL calculates the states

x = s0(1,:);
y = s0(2,:);
th = u(1:N);
A = u(N+1:2*N);
dt = u(2*N+1);

% velocities
Vvector = zeros(2,N);
Vguess = Vvector; Vguess(:,1) = [1;1];
Vriver = zeros(2,N);

% state propagation
for i = 1:N-1
    Vriver(1,i) = calcVriver(y(i));
    Vvector(:,i) = calcVvector(A(i),th(i),Vguess(:,i),Vriver(:,i));
    Vguess(:,i+1) = Vvector(:,i);
    x(i+1) = x(i) + Vvector(1,i)*dt;
    y(i+1) = y(i) + Vvector(2,i)*dt;
end
s = [x; y];



% % velocities
% Vvector = ones(2,N); Vvector(:,1) = [1,1];
% 
% % state propagation
% V = zeros(N,2); Vriver = zeros(N,1);
% for i = 1:N
%     Vriver(i) = calcVriver(y(i));
%     Vvector(:,i) = calcVvector(Vw,Vriver(i),A(i),th(i),Vvector(:,i));
%     x(i+1) = x(i) + Vvector(1,i)*dt;
%     y(i+1) = y(i) + Vvector(2,i)*dt;
% end
% s = [x; y];