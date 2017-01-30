function [  ] = plotStates( t, state, state_data, n_states, mean, sigma )
%PLOTSTATES Summary of this function goes here
%   t = time vector
%   state_date = time realizations of the state variable
%   n_states = number of state realizations to plot
%   mean = ensemble mean of the state variable for all time
%   sigm = ensemble standard dev of the state variable for all time

% f(state) = figure(state);
figure(state)
hold on
hx_1 = plot(t,state_data(:,1:n_states),'k');
hx_2 = plot(t,mean,'r',...
            t,mean+sigma,'b',...
            t,mean-sigma,'b',...
            'linewidth',3);
hold off
title('10 sample realizations and ensemble statistics (1000 realizations)')
statename = cellstr(['x-position';'x-velocity';'y-position']);
ylabel(statename(state))
xlabel('Time (s)')
legend(hx_2,'Calculated Mean','Calculated \sigma','location','southwest');

end

