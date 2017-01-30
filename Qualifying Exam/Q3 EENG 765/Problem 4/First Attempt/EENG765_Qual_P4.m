% Tim Coon 11/26/2014
% Q3 Qualifying Exam EENG765
% Problem #4
clear; close all; clc;

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
n_reals = 1000;         % number of realizations
t_span = 7;             % (sec) simulation time
t1 = 0:dt_x:t_span;     % dynamics time vector
n_states = length(t1);  % number of time steps in process
% shaping filters
sFilter_x = sqrt(2*sigma_x^2*beta_x);
sFilter_y = sqrt(2*sigma_yn^2*beta_y);
% State-Space matrices from hand calculations
fx = [0 1; 0 -beta_x];          % x-direction system matrix
fy = [0];                       % y-direction system matrix
nx_statevar = length(fx);       % number of x-direction state variables
ny_statevar = length(fy);       % number of y-direction state variables
n_statevar = nx_statevar + ny_statevar;   % total number of state variables
gx = [0; sFilter_x];            % x-direction input coeff matrix
gy = [sFilter_y];               % y-direction input coeff matrix
% calculate assumed-static parameters
F = get_F(fx, fy);                  % system matrix
G = get_G(gx, gy);                  % input coefficients
% [state transition matrix, discrete covariance matrix]
[phi,Qd] = get_phi_Qd(F,G,W,dt_x);  
% initial conditions
X = zeros(n_statevar,n_states,n_reals);

%% generate sample realizations of the dynamics
% generate filtered noise vectors for each realization
wd = zeros(n_statevar,n_states,n_reals);
for run = 1:n_reals
    wd(:,:,run) = discreteNoiseGenerator(Qd,n_sources,n_states,G);
end
for run = 1:n_reals
%     wd(:,:,run) = discreteNoiseGenerator(Qd, n_sources, n_states, G);
    for k = 2:n_states
        X(:,k,run) = phi*X(:,k-1,run) + wd(:,k-1,run);
    end
end

%% Plot real dynamics
% plot 10 realizations of each state, then plot the ensemble statistics of
% 10000 realizations of each state
for statevar = 1:2:n_statevar
    % use sqeeze() to isolate state data as a 2D matrix
    [s_mean(statevar,:), s_sigma(statevar,:)] = calcEnsembleStats(squeeze(X(statevar,:,:)));
    plotStates(t1,statevar,squeeze(X(statevar,:,:)),10,s_mean(statevar,:),s_sigma(statevar,:));
end

%% Measurement Realizations
% given measurement parameters
B_var = 2;                          % (m^2)   constant measurement bias
t2 = [0 3 5];                       % sample-time vector
n_meas = 3;                         % total number of measurements taken (predefined)
H = [0.4 0 2];                      % Measurement Matrix
n_sense = size(H,1);                % number of sensors and sensor noises
R = 0;                              % Measurement covariance matrix
% generate measurements from the first realization
X = squeeze(X(:,:,1));                  % use only the first realization
Z = zeros(n_sense,n_states);              % initialize measurement matrix
vd = R*randn(n_sense,n_states);           % persistent measurement noise/uncertainty (there is none)
B = sqrt(B_var)*randn;                  % constant bias on measuremnts
for i = 1:n_states
    Z(:,i) = H*X(:,i);
end

%% Estimation
xm0 = zeros(n_statevar,1);
Pm0 = zeros(n_statevar);            % start with Qd
[x_hat, x_std] = TC_KF(xm0,Pm0,phi,H,Qd,R,Z,t1);

%% Estimation Plots
figure()
suptitle({'Title';' '})
set(0,'Units','pixels')
sz = get(0,'ScreenSize');
set(gcf,'Position',[0 0 sz(3)/2 sz(4)])

% px ----------------------------------------------------------------------
subplot(311)
hold on
plot(t1,X(1,:),'b','linewidth',2)
plot(t1,x_hat(1,:),'r--','linewidth',2)
plot(t1,x_std(1,:),'k','linewidth',1)
plot(t1,-x_std(1,:),'k','linewidth',1)
hold off
ylabel('x-position (m)'); xlabel('time (s)');
hl1 = legend('px Truth','Estimate','St Dev','location','eastoutside');
pos_hl1 = get(hl1,'position'); pos_hl1(4) = 2*pos_hl1(4);
set(hl1,'position',pos_hl1)

% py ----------------------------------------------------------------------
subplot(312)
hold on
plot(t1,X(3,:),'b','linewidth',2)
plot(t1,x_hat(3,:),'r--','linewidth',2)
plot(t1,x_std(3,:),'k','linewidth',1)
plot(t1,-x_std(3,:),'k','linewidth',1)
hold off
ylabel('y-position (m)'); xlabel('time (s)');
hl2 = legend('py Truth','Estimate','St Dev','location','eastoutside');
pos_hl2 = get(hl2,'position'); pos_hl2(4) = 2*pos_hl2(4);
set(hl2,'position',pos_hl2)

% path --------------------------------------------------------------------
subplot(313)
hold on
plot(X(1,:),X(3,:),'b-*','linewidth',2)
plot(x_hat(1,:),x_hat(3,:),'r--*','linewidth',2)
hold off
ylabel('y-position (m)'); xlabel('x-position (m)');
hl3 = legend('True Path','Estimated Path','location','eastoutside');
pos_hl3 = get(hl3,'position'); pos_hl3(4) = 2*pos_hl3(4);
set(hl3,'position',pos_hl3)