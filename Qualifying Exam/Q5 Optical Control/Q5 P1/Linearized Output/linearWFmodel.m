function [  ] = linearWFmodel( X, beam )
%LINEARWFMODEL calculates linear WF from dynamics and plots results
%   X = perturbation state vectors (Eq. 113)
%   C_WF = linearized output matrix


%% generate output wavefront state vectors
wn = beam.C_WF * X;

% transform this into reference surface coordinate frame (use psi)
x_hat_rs = [1; 0; 0];
y_hat_rs = [0; 0; 1];
T_rs = [[x_hat_rs.'; y_hat_rs.'], zeros(2,4)];

%% Spot diagram
% nominal spot diagram
beam.plotNomSpots(beam);
% perturbed spot diagrams
for j = 1:length(X);
    for i = 1:beam.nrays
        gamma(:,i) = wn(4*i:4*i+3,j);
    end
    beam.plotSpots(beam,gamma);
end

%% OPD distribution plot


end

