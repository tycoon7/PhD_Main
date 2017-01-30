% Simple mirror system test
clear all; close all; clc;

% Mirror Parameters
f = .780;           % (m) focal length of OAP1
e = 1;              % (-) eccentricity of OAP1
psi = [0; 1; 0];    % (-) principal axis direction of OAP1
M = eye(3)-e^2*(psi*transpose(psi));  % surface dyadic
q = [0; 0; 0];      % mirror rotation point

if f ~= inf
    N_0 = -f*(1+e)*psi;
else
    N_0 = psi;
end

% Light path parameters
p = [-0.1; 1.5; 0];    % starting point of incident ray
i = [0; -1; 0];        % direction of incident ray

% find L1
L = findL(i,M,p,N_0);

% calc rho_1
rho = p+L*i;

% find direction of reflected ray
if f ~= inf
    N_vec = N_0 + M*rho;       % surface normal vector at intersection
    N_hat = -sign(dot(i,N_vec))*N_vec/norm(N_vec);
else
    N_vec = psi;
    N_hat = psi;
end
% N1_vec = N_01 + M1*rho_1;       % surface normal vector at intersection
% N1 = norm(N1_vec);              % magnitude of surface normal
% N1_hat = -sign(dot(i1,N1_vec))*N1_vec/N1;
R = eye(3)-2*N_hat*transpose(N_hat);

r = R*i;      % verified r1 as passing through focal point

% move this process into a function as much as possible for application to
% subsequent elements in the system

%% sensitivities
dr_di = drdi(N_hat,i,N_vec,M,L,R);
dr_da = drda(N_hat,i,N_vec,M);
dr_dt = drdt(e,i,r,N_hat,N_vec,M,p,q);
dr_dd = drdd(i,N_hat,N_vec,M);

dg_di = dgdi(L,R);
dg_da = dgda(R);
dg_dd = dgdd(r,i,N_hat);
dg_dt = dgdt(r,i,N_hat,p,q);

dL_di = zeros(3,1);
dL_da = zeros(3,1);
dL_dd = dLdd(r,i,N_hat);
dL_dt = dLdt(i,r,N_hat,p,q);

%% plots
figure(1)
hold on
plot(rho(1),rho(2),'g+','MarkerSize', 30)
plot(0,0,'.','MarkerSize',20)
plot(0,f,'+','MarkerSize',20)
quiver(p(1),p(2),i(1),i(2),'Color','r')
quiver(rho(1),rho(2),r(1),r(2),'Color','r')
hold off
xlabel('x'); ylabel('y');
axis equal