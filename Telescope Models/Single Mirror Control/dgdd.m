function [ dgdd ] = dgdd(r,i,N_hat)
%DGDD Summary of this function goes here
%   Detailed explanation goes here

P_r = eye(3) - r*r';
dgdd = P_r*i*N_hat'/dot(i,N_hat);   % eqn (72)

end

