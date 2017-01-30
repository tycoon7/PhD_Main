function [phi,Q_d] = get_phi_Qd(F,G,W,dt)
% Van Loan Method (Brown & Hwang sec 3.9)
% F = state-space system matrix
% G = state-space input matrix
% W = continuous noise-covariance matrix (contains PSDs)
%     (when using shaping filter notation of B&H sec 3.6, the white noise
%     scaling factors are in G, so W is the identity matrix because the
%     noise inputs are unity white noise.)
% dt = time step

NS = size(F,1);
NN = size(G,2);

A = [   -F           G*W*G';
     zeros(NS,NS)      F'  ] * dt;
B = expm(A);

B12 = B(1:NS,NS+1:2*NS);
B22 = B(NS+1:2*NS, NS+1:2*NS);

phi = B22';
Q_d = phi * B12;
    