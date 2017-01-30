function [ u_ZOH ] = make_ZOH_vector(u,nstates)
%MAKE_ZOH_VECTOR expands a vector of values by zero-order-hold
%   Take in a vector of values and extend each value to cover the entire
%   time step the value is held over. This was written for EENG765 qual
%   because the y-position input updates every second, but I want the
%   dynamics of the x-position to update faster.
% 
% u = control vector before expansion
% nstates = desired number of states

L_u = length(u);            % length of input vector
L_ZOH = nstates/L_u;        % length of ZOH step (1 sec)
u_ZOH = zeros(nstates,1);   % intitialize the output vector

for i = 1:L_u
    u_ZOH((i-1)*L_ZOH+1:i*L_ZOH) = u(i)*ones(L_ZOH,1);
end

% There are dimensional issues if nstates is not a multiple of L_u

end

