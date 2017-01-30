% CTEx Telescope
% close all;
clear; clc;

% Mirror Parameters
f1 = inf;                           % (m) focal length of SSM
f2 = 1.1633;                        % (m) focal length of OAP1
f3 = 0.1454;                        % (m) focal length of OAP2
f4 = inf;                           % (m) focal length of FSM
f5 = inf;                           % (m) focal length of beamsplitter  
f6 = inf;                           % (m) focal length of REF
e1 = 0;                             % (-) eccentricity of SSM
e2 = 1;                             % (-) eccentricity of OAP1
e3 = 1;                             % (-) eccentricity of OAP2
e4 = 0;                             % (-) eccentricity of FSM
e5 = 0;                             % (-) eccentricity of beamsplitter
e6 = 0;                             % (-) eccentricity of reference surface
V1 = [-1.20; 1.67; 0.00];           % (m) vertex, SSM
V2 = [-0.87; 0.45; 0.00];           % (m) vertex, OAP1
V3y = V2(2)+f2-f3;                  % calculate y-position of M3 vertex
V3 = [-0.87; V3y;  0.00];           % (m) vertex, OAP2
V4 = [-0.91; 0.20; 0.00];           % (m) vertex, FSM
V5 = [ 0.00; 0.20; 0.00];           % (m) vertex, BSp
V6 = [ 0.00; 0.00; 0.00];           % (m) vertex, RS
Thz1 = 225;                         % (deg)
Thz2 = 90;                          % (deg)
Thz3 = 90+.1;                          % (deg)
Thz4 = 45;                          % (deg)
Thz5 = 225;                         % (deg)
Thz6 = 90;                          % (deg)
psi1 = [cosd(Thz1); sind(Thz1); 0];        % (-) principal axis
psi2 = [cosd(Thz2); sind(Thz2); 0];        % (-) principal axis
psi3 = [cosd(Thz3); sind(Thz3); 0];        % (-) principal axis
psi4 = [cosd(Thz4); sind(Thz4); 0];        % (-) principal axis
psi5 = [cosd(Thz5); sind(Thz5); 0];        % (-) principal axis
psi6 = [cosd(Thz6); sind(Thz6); 0];        % (-) principal axis

e = [e1 e2 e3 e4 e5 e6];
f = [f1 f2 f3 f4 f5 f6];
V = [V1 V2 V3 V4 V5 V6];
psi = [psi1 psi2 psi3 psi4 psi5 psi6];
% Thz = [Thz1 Thz2 Thz3 Thz4 Thz5 Thz6];
Thz = [0 0 0 0 0 0];

P0 = [-1.50; 1.67; 0.00];
I0 = [1; 0; 0];

%% 

% build mirror objects
for s = 1:length(e)
    M(s) = MirrorClass(e(s),f(s),V(:,s),psi(:,s),rotz(-Thz(s)));
end

% build ray objects
nrays = 5;
for r = 1:nrays
    ray(r) = RayClass(r,P0,I0,M);
end

% calculate sensitivity matrices


%% plots

% overview
figure(1)
hold on
for s = 1:length(e)+1
    for r = 1:nrays
        % intersection points
        plot3(ray(r).seg(s).P(1),ray(r).seg(s).P(2),ray(r).seg(s).P(3),'g+','MarkerSize', 30);
        % ray direction
        if s ~= length(e)+1
%         quiver3(ray(r).seg(s).P(1),ray(r).seg(s).P(2),ray(r).seg(s).P(3),...
%                 ray(r).seg(s).I(1),ray(r).seg(s).I(2),ray(r).seg(s).I(3),'Color','r')
        arrow([ray(r).seg(s).P(1),ray(r).seg(s).P(2),ray(r).seg(s).P(3)],...
              [ray(r).seg(s+1).P(1),ray(r).seg(s+1).P(2),ray(r).seg(s+1).P(3)],...
              'EdgeColor','r','FaceColor','r')
        end
    end
end
for j=1:length(e)
    plot3(V(1,j),V(2,j),V(3,j),'.','MarkerSize',20);
end
hold off
axis equal
xlabel('x'); ylabel('y');
% view(3)


% reference surface spot diagram
figure(2)
hold on
for r = 1:nrays
    plot(ray(r).seg(end).P(1),ray(r).seg(end).P(3),'.','MarkerSize',20);
end
hold off
ylim([-.05 .05]);
xlim([-.05 .05]);
