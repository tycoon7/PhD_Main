function [ spotSize ] = SingleMirror_SpotSize( alpha, V1)
%SINGLEMIRROR_SPOTSIZE returns spot size for given theta misalignment and
%vertex location
%   alpha = rotation of primary mirror
%   V1 = position of the vertex, initially located at the origin
% In the simulation code, the vertex location is found wrt the mirror
% centroid. This raytrace code assumes the vertex is at the origin so the
% values input are deviations from the at rest position.

global f_par D Pc

% Mirror Parameters
dR = rotz(rad2deg(alpha));
f1 = f_par;                         % (m) focal length of parabolic mirror
f2 = inf;                           % (m) focal length of reference surface
e1 = 1;                             % (-) eccentricity of par
e2 = 0;                             % (-) eccentricity of rs
V1 = [V1(1); V1(2);  0.00];         % (m) vertex, par
V2 = [0.00;  f1;  0.00];            % (m) vertex, rs
q1 = dR*[Pc(1); Pc(2); 0.00];       % (m) rotation point, par
q2 = V2;                            % (m) rotation point, rs
% Thz1 = 0;                           % (deg)
% Thz2 = 0;                           % (deg)
psi1 = [0;  1; 0];                  % (-) principal axis of par
psi2 = [0; -1; 0];                  % (-) principal axis of rs

% Assemble parameters
e = [e1 e2];
f = [f1 f2];
V = [V1 V2];
q = [q1 q2];
psi = [dR*psi1 psi2];
% Thz = [Thz1 Thz2];

P0 = [-0.25; 1.5; 0];
I0 = [0; -1; 0];

%% 
% build mirror objects
for s = 1:length(e)
    M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),eye(3));
%     M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),rotz(-Thz(s)));
end

% build ray objects
% nrays = 3;
% for r = 1:nrays
%     ray(r) = RayClass(r,P0,I0,D,M);
% end

% build beam objects
beam = BeamClass(D,5,P0,I0,M);

%% output matrix
% This is if we need the sensitivities. Functions are with the original
% code for the single mirror control

%% maximum ray intersection distance from chief ray on the detector plane
d = zeros(1,length(e));
for j = 2:length(e)
    d(j) = norm(beam.ray(j).seg(end).P(1)-beam.ray(1).seg(end).P(1),beam.ray(j).seg(end).P(3)-beam.ray(1).seg(end).P(3));
end
spotSize = max(d);

%% plots

% overview
figure(1)
hold on
parfor s = 1:length(e)-1
    for r = 1:beam.nrays
        % intersection points
        plot3(beam.ray(r).seg(s).P(1),beam.ray(r).seg(s).P(2),beam.ray(r).seg(s).P(3),'g+','MarkerSize', 30);
        % ray direction
        if s ~= length(e)+1
        quiver3(beam.ray(r).seg(s).P(1),beam.ray(r).seg(s).P(2),beam.ray(r).seg(s).P(3),...
                beam.ray(r).seg(s).I(1),beam.ray(r).seg(s).I(2),beam.ray(r).seg(s).I(3),'Color','r')
        end
    end
end
for j=1:length(e)
    plot3(V(1,j),V(2,j),V(3,j),'.','MarkerSize',20);
end
hold off
axis equal
xlabel('x'); ylabel('y');
view(3)


% reference surface spot diagram
figure(2)
hold on
for r = 1:beam.nrays
    plot(beam.ray(r).seg(end).P(1),beam.ray(r).seg(end).P(3),'.','MarkerSize',20);
end
hold off
xlabel('x'); ylabel('z');
ylim([-.05 .05]);
xlim([-.05 .05]);


