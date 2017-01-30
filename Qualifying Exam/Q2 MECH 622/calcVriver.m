function [ Vriver ] = calcVriver( y )
%CALCVRIVER calculates the velocity of the river
%   Detailed explanation goes here

global Vmax w

% Vriver = 0;
% if y < 0 || y > 200
%     Vriver = 10;
% else
%     Vriver = 4*Vmax*y*(w-y)/w^2;

Vriver = abs(4*Vmax*y*(w-y)/w^2);

end

