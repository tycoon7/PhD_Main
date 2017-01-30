function [ drdt ] = drdt(e,i,r,N_hat,N_vec,M,p,q)
%DRDT calculate sensitivity drdt 
%   Detailed explanation goes here

P_N = -skew(N_hat)*skew(N_hat);

% for flat mirrors, the equation is simplified. Check e
if e == 0
    drdt = -2*skew(r)*P_N;              % eqn (70)
else
    drdN = -2*(dot(N_hat,i)*eye(3)+N_hat*i');   % eqn (43)
    dNdt = -skew(N_hat);                        % eqn (66)
    dNdrho = -sign(dot(i,N_vec))*(1/norm(N_vec))*P_N*M;    % eqn (21)
    drhodd = -(eye(3)-i*N_hat'/dot(i,N_hat));   % eqn (63)
    drdd = drdN*dNdrho*drhodd;                  % eqn (64)
    dddt = -skew(p-q);                          % eqn (68)
    drdt = drdN*dNdt + drdd*dddt;               % eqn (69)
end

