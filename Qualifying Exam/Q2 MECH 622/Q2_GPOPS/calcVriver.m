function [ Vriver ] = calcVriver( y )
%CALCVRIVER calculates the x-velocity of the river
%   Detailed explanation goes here

global Vmax w

% Vriver = 0;
Vriver = 4*Vmax*y*(w-y)/w^2;

end

