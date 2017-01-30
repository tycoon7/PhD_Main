function [ drdd ] = drdd(i,N_hat,N_vec,M)
%DRDD calculate sensitivity drdd
%   Detailed explanation goes here

P_N = -skew(N_hat)*skew(N_hat);             % eqn (70)
drdN = -2*(dot(N_hat,i)*eye(3)+N_hat*i');   % eqn (43)
dNdrho = -sign(dot(i,N_vec))*(1/norm(N_vec))*P_N*M;    % eqn (21)
drhodd = -(eye(3)-i*N_hat'/dot(i,N_hat));   % eqn (63)
drdd = drdN*dNdrho*drhodd;                  % eqn (64)

end

