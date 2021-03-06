function [ drdi ] = drdi( N_hat, i, N_vec, M, L, R )
%DRDI calculate sensitivity, drdi
%   Detailed explanation goes here

drdN = -2*(dot(N_hat,i)*eye(3)+N_hat*i');  % eqn (43)
P_N = -skew(N_hat)*skew(N_hat);
dNdrho = -sign(dot(i,N_vec))*(1/norm(N_vec))*P_N*M;    % eqn (21)
drhodi = L*(eye(3)-i*N_hat'/dot(i,N_vec)); % eqn (51)
drdi = R + drdN*dNdrho*drhodi;   % eqn (52)
end

