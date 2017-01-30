% Simple mirror control system test
% This system composes a single parabolic mirror and a reference surface
% with a series of axial rays at the input. 
% close all;
clear all; clc;

% Mirror Parameters
f1 = .780;                          % (m) focal length of parabolic mirror
f2 = inf;                           % (m) focal length of reference surface
e1 = 1;                             % (-) eccentricity of par
e2 = 0;                             % (-) eccentricity of rs
V1 = [0.00; -0.78; 0.00];           % (m) vertex, par
V2 = [0.00;  0.00; 0.00];           % (m) vertex, rs
Thz1 = 0;                           % (deg)
Thz2 = 0;                           % (deg)
psi1 = [0;  1; 0];                  % (-) principal axis of par
psi2 = [0; -1; 0];                  % (-) principal axis of rs
dR = rotz(0.1);

e = [e1 e2];
f = [f1 f2];
V = [V1 V2];
psi = [dR*psi1 psi2];
Thz = [Thz1 Thz2];

P0 = [-0.25; 1.5; 0];
I0 = [0; -1; 0];

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
% dr_di = drdi(N_hat,i,N_vec,M,L,R);
% dr_da = drda(N_hat,i,N_vec,M);
% dr_dt = drdt(e,i,r,N_hat,N_vec,M,p,q);
% dr_dd = drdd(i,N_hat,N_vec,M);
% 
% dg_di = dgdi(L,R);
% dg_da = dgda(R);
% dg_dd = dgdd(r,i,N_hat);
% dg_dt = dgdt(r,i,N_hat,p,q);
% 
% dL_di = zeros(3,1);
% dL_da = zeros(3,1);
% dL_dd = dLdd(r,i,N_hat);
% dL_dt = dLdt(i,r,N_hat,p,q);

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
        quiver3(ray(r).seg(s).P(1),ray(r).seg(s).P(2),ray(r).seg(s).P(3),...
                ray(r).seg(s).I(1),ray(r).seg(s).I(2),ray(r).seg(s).I(3),'Color','r')
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