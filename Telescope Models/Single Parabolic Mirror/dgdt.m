function [ dgdt ] = dgdt(r,i,N_hat,p,q)
%DGDT Calculate Sensitivite dgdt
%   Detailed explanation goes here

P_r = eye(3) - r*r';
dgdd = P_r*i*N_hat'/dot(i,N_hat);   % eqn (72)
dddt = -skew(p-q);                  % eqn (68)
dgdt = dgdd*dddt;

end

