function [ C_WF ] = linOutputMatrix()
%LINOUTPUTMATRIX returns the output matrix to approximate the wavefront
%   C_WF = WaveFront model linear output matrix

global f_par Pc

% Mirror Parameters
f1 = f_par;                         % (m) focal length of parabolic mirror
f2 = inf;                           % (m) focal length of reference surface
e1 = 1;                             % (-) eccentricity of par
e2 = 0;                             % (-) eccentricity of rs
V1 = [0.00; 0.00;  0.00];         % (m) vertex, par
V2 = [0.00;  f1;  0.00];            % (m) vertex, rs
q1 = dR*[Pc(1); Pc(2); 0.00];       % (m) rotation point, par
q2 = V2;                            % (m) rotation point, rs
psi1 = [0;  1; 0];                  % (-) principal axis of par
psi2 = [0; -1; 0];                  % (-) principal axis of rs

% Assemble parameters
e = [e1 e2];
f = [f1 f2];
V = [V1 V2];
q = [q1 q2];
psi = [psi1 psi2];

P0_chief = [-0.25; 1.5; 0];
I0 = [0; -1; 0];

%% 
% build mirror objects
for s = 1:length(e)
    M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),eye(3));
%     M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),rotz(-Thz(s)));
end

% build beam objects
rd = 3;         % ray density
beam = BeamClass(D,rd,P0_chief,I0,M);

