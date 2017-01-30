function [ x, x_std ] = TC_KF_P4( xm0, Pm0, phi, H, Q, R, n_states, t )
%TC_KF generates state vector estimate at all states
%   Inputs:
%   xm0 = (nx1) initial state vector minus estimate
%   Pm0 = (nxn) initial error covariance matrix minus estimate
%   phi = (nxn) state transition matrix
%   H   = (mxn) output matrix
%   Q   = (nxn) process noise covariance matrix
%   R   = (mxm) measurement noise covariance matrix
%   n_states   = number of states
%   t   = (1xk) time vector
% 
%   Outputs:
%   x_hat = (n,k) columns are updated state vectors at each state
%   x_std = (n,k) updated standard deviations
%   
%   Takes in process parameters (x_minus0, P_minus0, phi, H, Q, R) assumed
%   to be constant.
%   
%   n_meas in not necessarily the same as the n_states of the simulation,

n_var = length(xm0);
lP = length(Pm0);
Pm = zeros(lP,lP,n_states);
Pm(:,:,1) = Pm0;
xm = zeros(n_var,n_states);
xm(:,1) = xm0;
if isempty(R)
    R = 0;
end

for k = 1:n_states
    if t(k)==3    % measurement times
        Z = 3.5;
        % step 1 - Kalman Weight
%         K = (Pm(:,:,k)*H')/(H*Pm(:,:,k)*H'+R);          % (4.2.17)
        K = (Pm(:,:,k)*H')*inv((H*Pm(:,:,k)*H'+R));          % (4.2.17)
        % step 2 - update estimate
        x(:,k) = xm(:,k)+K*(Z-H*xm(:,k));          % (4.2.8)
        % step 3 - update error covariance
        P(:,:,k) = (eye(lP)-K*H)*Pm(:,:,k);             % (4.2.22)
        % step 4 - propagate forward for initial estimate of next state
        xm(:,k+1) = phi*x(:,k);                         % (4.2.23)
        Pm(:,:,k+1) = phi*P(:,:,k)*phi' + Q;            % (4.4.25)
        
    elseif t(k) == 5
        Z = 5;
        % step 1 - Kalman Weight
%         K = (Pm(:,:,k)*H')/(H*Pm(:,:,k)*H'+R);          % (4.2.17)
        K = (Pm(:,:,k)*H')*inv((H*Pm(:,:,k)*H'+R));          % (4.2.17)
        % step 2 - update estimate
        x(:,k) = xm(:,k)+K*(Z-H*xm(:,k));          % (4.2.8)
        % step 3 - update error covariance
        P(:,:,k) = (eye(lP)-K*H)*Pm(:,:,k);             % (4.2.22)
        % step 4 - propagate forward for initial estimate of next state
        xm(:,k+1) = phi*x(:,k);                         % (4.2.23)
        Pm(:,:,k+1) = phi*P(:,:,k)*phi' + Q;            % (4.4.25)
        
    else
        % project estimates forward
        x(:,k) = xm(:,k);
        P(:,:,k) = Pm(:,:,k);
        % step 4 - propagate forward for initial estimate of next state
        xm(:,k+1) = phi*x(:,k);                         % (4.2.23)
        Pm(:,:,k+1) = phi*P(:,:,k)*phi' + Q;            % (4.4.25)
    end
    x_std(:,k) = sqrt(diag(P(:,:,k)));
end

end

