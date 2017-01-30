function [ J ] = Cost_MECH622_Qual( u, N, s0  )
%COST_MECH622_QUAL calculate value of the cost function
%   The last entry in the control vector is the time step. The cost
%   function is simply the value of time final

tf = u(end)*N;

J = tf;

end

