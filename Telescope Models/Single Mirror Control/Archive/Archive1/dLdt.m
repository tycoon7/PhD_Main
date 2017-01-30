function [ dLdt ] = dLdt(i,r,N_hat,p,q)
%DLDT Calculate Sensitivity dLdt
%   Detailed explanation goes here

dLdd = -((1-dot(r,i))/dot(i,N_hat))*N_hat';     % eqn (78)
dddt = -skew(p-q);                              % eqn (68)
dLdt = dLdd*dddt;                               % eqn (79)

end

