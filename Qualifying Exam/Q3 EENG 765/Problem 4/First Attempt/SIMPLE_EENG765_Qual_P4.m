% Tim Coon 11/26/2014
% Q3 Qualifying Exam EENG765
% Problem #4
clear all; close all; clc;

%% Dynamics Realizations
% given parameters
dt_x = 1;               % (sec) propagation time = time step duration
dt_y = 1;               
beta_x = 1/dt_x;        % inverse of the time constant
beta_y = 1/dt_y;              
sigma_x = 1;            % (m)
sigma_yn = 3;           % (m)
n_sources = 2;          % number of unity white noise sources
W = eye(n_sources);     % PSD matrix for Van Loan Method (identity for UWN)
t_span = 7;             % (sec) simulation time
t1 = 0:dt_x:t_span;     % dynamics time vector
n_states = length(t1);  % number of time steps in process
% shaping filters
sFilter_x = sqrt(2*sigma_x^2*beta_x);
sFilter_y = sqrt(2*sigma_yn^2*beta_y);
% State-Space matrices from hand calculations
fx = [0 1; 0 -beta_x];          % x-direction system matrix
fy = [0];                       % y-direction system matrix
n_statevar = length(fx) + length(fy);   % total number of state variables
gx = [0; sFilter_x];            % x-direction input coeff matrix
gy = [sFilter_y];               % y-direction input coeff matrix
% calculate assumed-static parameters
F = get_F(fx, fy);                  % system matrix
G = get_G(gx, gy);                  % input coefficients
% [state transition matrix, discrete covariance matrix]
[phi,Qd] = get_phi_Qd(F,G,W,dt_x);  

%% Measurement Realizations
H = [0.4 0 2];                      % Measurement Matrix
n_sense = size(H,1);                % number of sensors and sensor noises
R = 0;                              % Measurement covariance matrix

%% Estimation
xm0 = zeros(n_statevar,1);
Pm0 = zeros(n_statevar);            % start with Qd
[x_hat, x_std] = TC_KF_P4(xm0,Pm0,phi,H,Qd,R,n_states,t1);

%% Estimation Plots
figure()
suptitle({'Estimate of Ant Location';' '})
set(0,'Units','pixels')
sz = get(0,'ScreenSize');
set(gcf,'Position',[0 0 sz(3)/2 sz(4)])
% px ----------------------------------------------------------------------
subplot(311)
hold on
plot(t1,x_hat(1,:),'r--','linewidth',2)
plot(t1,x_std(1,:),'k','linewidth',1)
plot(t1,-x_std(1,:),'k','linewidth',1)
hold off
ylabel('x-position (m)'); xlabel('time (s)');
hl1 = legend('Estimate','St Dev','location','eastoutside');
pos_hl1 = get(hl1,'position'); pos_hl1(4) = 2*pos_hl1(4);
set(hl1,'position',pos_hl1)
% py ----------------------------------------------------------------------
subplot(312)
hold on
plot(t1,x_hat(3,:),'r--','linewidth',2)
plot(t1,x_std(3,:),'k','linewidth',1)
plot(t1,-x_std(3,:),'k','linewidth',1)
hold off
ylabel('y-position (m)'); xlabel('time (s)');
hl2 = legend('Estimate','St Dev','location','eastoutside');
pos_hl2 = get(hl2,'position'); pos_hl2(4) = 2*pos_hl2(4);
set(hl2,'position',pos_hl2)
% path --------------------------------------------------------------------
subplot(313)
hold on
plot(x_hat(1,:),x_hat(3,:),'r--*','linewidth',2)
hold off
ylabel('y-position (m)'); xlabel('x-position (m)');
hl3 = legend('Estimated Path','location','eastoutside');
pos_hl3 = get(hl3,'position'); pos_hl3(4) = 2*pos_hl3(4);
set(hl3,'position',pos_hl3)
