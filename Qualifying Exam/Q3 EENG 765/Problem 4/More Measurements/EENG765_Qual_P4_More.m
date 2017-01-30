% Tim Coon 11/26/2014
% Q3 Qualifying Exam EENG765
% Problem #4
clear; close all; clc;
set(0,'defaulttextinterpreter','latex')

%% Dynamics Realizations
% given parameters
dt_x = 0.1;               % (sec) propagation time = time step duration
dt_y = 1;               
beta_x = 1/dt_x;        % inverse of the time constant
beta_y = 1/dt_y;              
sigma_x = 1;            % (m)
sigma_yn = 3;           % (m)
% n_sources = 4;          % number of unity white noise sources
% W = eye(n_sources);     % PSD matrix for Van Loan Method (identity for UWN)
n_reals = 1000;         % number of realizations
t_span = 7;             % (sec) simulation time
t1 = 0:dt_x:t_span;     % dynamics time vector
n_states = length(t1);  % number of time steps in process
% given measurement parameters
b_var = 2;                  % (m^2) constant measurement bias variation
t2 = [3 5];               % sample-time vector
n_meas = length(t2);        % total number of measurements taken (given)
H = [0.4 0 2 1];            % Measurement Matrix
n_sense = size(H,1);        % number of sensors and sensor noises
R = [];                  % Measurement covariance matrix

%% Calculated values
% shaping filters
sFilter_x = sqrt(2*sigma_x^2*beta_x);
sFilter_y = sqrt(2*sigma_yn^2*beta_y);
% State-Space matrices from hand calculations
fx = [0 1; 0 -beta_x];          % x-direction system matrix
fy = [0];                       % y-direction system matrix
fb = [0];                       % bias system matrix
gx = [0 0; sFilter_x 0];        % x-direction input coeff matrix
gy = [0 sFilter_y];             %#ok<*NBRAK> % y-direction input coeff matrix
gb = [0 0];                     % bias input
% calculate assumed-static parameters
F = matrix_concat(fx, fy, fb);  % system matrix
n_statevar = size(F,1);
G = [gx; gy; gb];  % input coefficients, ref Maybeck (5-123)
n_sources = size(G,2);
W = eye(size(G,2));    % PSD matrix for Van Loan Method (identity for UWN)
% [state transition matrix, discrete covariance matrix]
[phi,Qd] = get_phi_Qd(F,G,W,dt_x);  
% initial conditions
X0 = zeros(n_statevar,1,n_reals);  % initial states
b0 = sqrt(b_var)*randn;            % constant bias on measuremnts
X0(4,1,:) = b0;               % initial value of bias for all realizations

%% generate sample realizations of the dynamics
% generate filtered noise vectors for each realization
wd = zeros(n_statevar,n_states,n_reals);
for run = 1:n_reals
    wd(:,:,run) = sqrt(Qd)*[zeros(1,n_states); randn(2,n_states); zeros(1,n_states)];
%     wd(:,:,run) = discreteNoiseGenerator(Qd,n_sources,n_states,G);
end
X = zeros(n_statevar,n_states,n_reals);  % initialize state realizations matrix
X(:,1,:) = X0;
for run = 1:n_reals
%     wd(:,:,run) = discreteNoiseGenerator(Qd, n_sources, n_states, G);
    for k = 2:n_states
        X(:,k,run) = phi*X(:,k-1,run) + wd(:,k-1,run);
    end
end

%% Plot real dynamics and their statistics
% plot 10 realizations of each state, then plot the ensemble statistics of
% 1000 realizations of each state
for statevar = 1:2:n_statevar
    % use sqeeze() to isolate state data as a 2D matrix
    [s_mean(statevar,:), s_sigma(statevar,:)] = calcEnsembleStats(squeeze(X(statevar,:,:)));
    plotStates(t1,statevar,squeeze(X(statevar,:,:)),10,s_mean(statevar,:),s_sigma(statevar,:));
end

%% Measurement Realizations
% generate measurements from the first realization
X = squeeze(X(:,:,1));                  % use only the first realization
Z = zeros(n_sense,n_states);            % initialize measurement matrix
% vd = R*randn(n_sense,n_states);         % persistent measurement noise/uncertainty (there is none)
for i = 1:n_states
    Z(:,i) = H*X(:,i);
end

%% Estimation
xm0 = zeros(n_statevar,1);
Pm0 = zeros(n_statevar);       % start known exactly
Pm0(end,end) = b_var;
[x_hat, x_std] = TC_KF(xm0,Pm0,phi,H,Qd,R,Z,t1);

%% Estimation Plots
figure()
suptitle({'Estimate of Ant Location';' '})
set(0,'Units','pixels')
sz = get(0,'ScreenSize');
set(gcf,'Position',[0 0 sz(3)/2 sz(4)])
% px ----------------------------------------------------------------------
subplot(411)
hold on
plot(t1,X(1,:),'b-','linewidth',2)
plot(t1,x_hat(1,:),'r--','linewidth',2)
stairs(t1,x_hat(1,:)+x_std(1,:),'k','linewidth',1)
stairs(t1,x_hat(1,:)-x_std(1,:),'k','linewidth',1)
hold off
ylabel('x-position (m)'); xlabel('time (s)');
legend('Estimate','St Dev','location','eastoutside');
% vx ----------------------------------------------------------------------
subplot(412)
hold on
plot(t1,X(2,:),'b-','linewidth',2)
plot(t1,x_hat(2,:),'r--','linewidth',2)
stairs(t1,x_hat(2,:)+x_std(2,:),'k','linewidth',1)
stairs(t1,x_hat(2,:)-x_std(2,:),'k','linewidth',1)
hold off
ylabel('x-velocity (m)'); xlabel('time (s)');
legend('Estimate','St Dev','location','eastoutside');
% py ----------------------------------------------------------------------
subplot(413)
hold on
plot(t1,X(3,:),'b-','linewidth',2)
plot(t1,x_hat(3,:),'r--','linewidth',2)
stairs(t1,x_hat(3,:)+x_std(3,:),'k','linewidth',1)
stairs(t1,x_hat(3,:)-x_std(3,:),'k','linewidth',1)
hold off
ylabel('y-position (m)'); xlabel('time (s)');
legend('Estimate','St Dev','location','eastoutside');
% b  ----------------------------------------------------------------------
subplot(414)
hold on
plot(t1,x_hat(4,:),'r--','linewidth',2)
stairs(t1,x_hat(4,:)+x_std(4,:),'k','linewidth',1)
stairs(t1,x_hat(4,:)-x_std(4,:),'k','linewidth',1)
hold off
ylabel('Sensor Bias'); xlabel('time (s)');
legend('Estimate','St Dev','location','eastoutside');

% path --------------------------------------------------------------------
figure()
subplot(3,1,1:2)
hold on
plot(x_hat(1,:),x_hat(3,:),'b-*','linewidth',2)
plot(X(1,:),X(3,:),'r--o','linewidth',2)
hold off
ylabel('y-position (m)'); xlabel('x-position (m)');
title({'Paths'; ''});
legend('Estimate','Actual')
% axis equal
% xlim([-1 1.5]); ylim([0 2.5]);
ax = subplot(3,1,3);
set(ax,'visible','off')
x_h = num2str(x_hat,'%10.4f');
text(0.2,0.3,x_h)
text(0.1,0.3,'$\hat{x} = $')