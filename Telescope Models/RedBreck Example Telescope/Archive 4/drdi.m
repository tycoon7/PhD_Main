function [ drdi ] = drdi( ray, seg, mir )
%DRDI calculate sensitivity, drdi
%   Detailed explanation goes here

drdN = -2*(dot(ray.seg.N_hat,ray.seg.i)*eye(3)+ray.seg.N_hat*ray.seg.i');  % eqn (43)
P_N = -skew(ray.seg.N_hat)*skew(ray.seg.N_hat);
dNdrho = -sign(dot(ray.seg.i,ra.seg.N_vec))*(1/ray.seg.N_vec)*P_N*mir.M;    % eqn (21)
drhodi = ray.seg.L*(eye(3)-ray.seg.i*ray.seg.N_hat'/dot(ray.seg.i,ra.seg.N_vec)); % eqn (51)
drdi = seg.R + drdN*dNdrho*drhodi;   % eqn (52)
end

