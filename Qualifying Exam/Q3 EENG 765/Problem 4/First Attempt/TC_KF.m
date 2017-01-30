function [ x, x_std ] = TC_KF( xm0, Pm0, phi, H, Q, R, z, t )
%TC_KF generates state vector estimate at all states
%   Inputs:
%   xm0 = (nx1) initial state vector minus estimate
%   Pm0 = (nxn) initial error covariance matrix minus estimate
%   phi = (nxn) state transition matrix
%   H   = (mxn) output matrix
%   Q   = (nxn) process noise covariance matrix
%   R   = (mxm) measurement noise covariance matrix
%   z   = (mxk) measurement vector
%   k   = number of states
%   t   = (1xk) time vector
%   Outputs:
%   x_hat = (n,k) columns are updated state vectors at each state
%   x_std = (n,k) updated standard deviations

n_var = length(xm0);
[~, n_states] = size(z);
lP = length(Pm0);
Pm = zeros(lP,lP,n_states);
Pm(:,:,1) = Pm0;
xm = zeros(n_var,n_states);
xm(:,1) = xm0;

for k = 1:n_states
    if t(k)==3 || t(k)==5    % measurement times
        % step 1 - Kalman Weight
%         K = (Pm(:,:,k)*H')/(H*Pm(:,:,k)*H'+R);          % (4.2.17)
        K = (Pm(:,:,k)*H')*inv((H*Pm(:,:,k)*H'+R));          % (4.2.17)
        % step 2 - update estimate
        x(:,k) = xm(:,k)+K*(z(:,k)-H*xm(:,k));          % (4.2.8)
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

