function [ dLdd ] = dLdd(r,i,N_hat)
%DLDD Calculate Sensitivity dLdd
%   Detailed explanation goes here

dLdd = -((1-dot(r,i))/dot(i,N_hat))*N_hat';     % eqn (78)

end

