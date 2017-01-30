SingleMirror_SpotSize_Lin( alpha, V1)
%SINGLEMIRROR_SPOTSIZE returns spot size for given theta misalignment and
%vertex location
%   alpha = rotation of primary mirror
%   V1 = position of the vertex, initially located at the origin
% In the simulation code, the vertex location is found wrt the mirror
% centroid. This raytrace code assumes the vertex is at the origin so the
% values input are deviations from the at rest position.

% global f_par D Pc
alpha = 0;
V1 = [0;0]
f_par = 1;
D = 0.2;
Pc = [-0.3; 0.1];

% Mirror Parameters
dR = rotz(rad2deg(alpha));
f1 = f_par;                         % (m) focal length of parabolic mirror
f2 = inf;                           % (m) focal length of reference surface
e1 = 1;                             % (-) eccentricity of par
e2 = 0;                             % (-) eccentricity of rs
V1 = [V1(1); V1(2);  0.00];         % (m) vertex, par
V2 = [0.00;  f1;  0.00];            % (m) vertex, rs
% q1 = dR*[Pc(1); Pc(2); 0.00];       % (m) rotation point, par
q1 = V1;
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

P0_chief = [-0.25; 1.5; 0];
I0 = [0; -1; 0];

%% 
% build mirror objects
for s = 1:length(e)
    M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),eye(3));
%     M(s) = MirrorClass(e(s),f(s),V(:,s),q(s),psi(:,s),rotz(-Thz(s)));
end

% % build ray objects
%     % build beam: find the coordinates circular grid start points
%     R = D/2;
%     rd = 3;            % ray density
%     [r,phi] = meshgrid(linspace(R/rd,R,rd),linspace(2*pi/rd,2*pi,rd));
%     X = r.*sin(phi);
%     Y = zeros(size(r));
%     Z = r.*cos(phi);
%     origin = [0; 0; 0];
%     rayStartGrid = [origin, [X(:)';Y(:)';Z(:)']];
% P0 = bsxfun(@plus,P0_chief,rayStartGrid);
% nrays = length(rayStartGrid);
% for nr = 1:nrays
%     ray(nr) = RayClass(nr,P0(:,nr),I0,M);
% end

% build beam objects
rd = 5;         % ray density
beam = BeamClass(D,rd,P0_chief,I0,M);

%% output matrix
% This is if we need the sensitivities. Functions are with the original
% code for the single mirror control

%% maximum ray intersection distance from chief ray on the detector plane
% d = zeros(1,length(e));
% for j = 2:length(e)
%     d(j) = norm(ray(nr).seg(end).P(1)-ray(1).seg(end).P(1),ray(nr).seg(end).P(3)-ray(1).seg(end).P(3));
% end
% spotSize = max(d);

%% plots

% overview
figure(1)
hold on
for s = 1:length(e)+1
    for nr = 1:beam.nrays
        % intersection points
        plot3(beam.ray(nr).seg(s).P(1),beam.ray(nr).seg(s).P(2),beam.ray(nr).seg(s).P(3),'g+','MarkerSize', 30);
        % ray direction
        if s ~= length(e)+1
        quiver3(beam.ray(nr).seg(s).P(1),beam.ray(nr).seg(s).P(2),beam.ray(nr).seg(s).P(3),...
                beam.ray(nr).seg(s).I(1),beam.ray(nr).seg(s).I(2),beam.ray(nr).seg(s).I(3),'Color','r')
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
for nr = 1:beam.nrays
    plot(beam.ray(nr).seg(end).P(1),beam.ray(nr).seg(end).P(3),'.','MarkerSize',20);
end
hold off
xlabel('x'); ylabel('z');
% ylim([-.05 .05]);
% xlim([-.05 .05]);


