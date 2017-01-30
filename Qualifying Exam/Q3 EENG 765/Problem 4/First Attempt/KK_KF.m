function [ x_hat_out, x_std_out ] = KK_KF( x0, P0, Phi_d, Q_d, H, R, z_meas, t1, t2 )
%KK_KF Summary of this function goes here
%   Detailed explanation goes here

[ti, ia, ib] = intersect(t1, t2);    % Locate the time indices where we want to do a measurement

n_states = length(x0);
NT = size(t1,1);
x_hat_minus = zeros(n_states,NT);
x_hat_plus = zeros(n_states,NT);
P_minus = zeros(n_states,n_states,NT);
P_plus = zeros(n_states,n_states,NT);

x_hat_minus(:,:,1) = x0;
x_hat_plus(:,:,1) = x0;
P_minus(:,:,1) = P0;
P_plus(:,:,1) = P0;

x_hat_out = zeros(n_states,2*NT);
x_std_out = zeros(n_states,2*NT);
x_hat_out(1:n_states,1) = x0;
x_std_out(1:n_states,1:2) = sqrt(P0);
% t_out = zeros(2*NT,1);
% t_out(1:2:2*NT) = t;
% t_out(2:2:2*NT) = t;

for ii=2:NT
    % Propagate state estimate and covariance from time ii to ii+1
    x_hat_minus(ii) = Phi_d * x_hat_plus(ii-1);
    P_minus(ii) = Phi_d * P_plus(ii-1) * Phi_d' + Q_d;
    
    x_hat_out(2*ii-1) = x_hat_minus(ii);
    x_std_out(2*ii-1) = sqrt(P_minus(ii));
    
    % If there is an update, incorporate it.  Else just transfer the state
    % esimtate
    idx_meas = find(ii == ia);
    if idx_meas
        % Update
        residual = z_meas(idx_meas) - H*x_hat_minus(ii);
        res_cov = H*P_minus(ii)*H' + R;
        K = P_minus(ii)*H'*inv(res_cov);
        
        x_hat_plus(ii) = x_hat_minus(ii) + K*residual;
        P_plus(ii) = P_minus(ii) - K*H*P_minus(ii);
    else
        x_hat_plus(ii) = x_hat_minus(ii);
        P_plus(ii) = P_minus(ii);
    end
    
    x_hat_out(2*ii) = x_hat_plus(ii);
    x_std_out(2*ii) = sqrt(P_plus(ii));
end
end

