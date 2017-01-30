% This file passes the angle of the mirrorDynamics function to the
% SingleMirror_SpotSize function to calculate the spot size for the
% duration of the simulation. The spot size is plotted for different focal
% lengths.
clear; close all; clc;

% load('Example.mat')
% openfig('tracking.fig')

f = [0.70 0.80];

[t,thetaZ,u,w] = mirrorDynamics();
% 
for c1 = 1:length(f)
    for c2 = 1:length(t)
        spotSize(c2,c1) = SingleMirror_SpotSize(thetaZ(c2),f(c1));
%         cost(c2,c1) = 1e5*spotSize(c2,c1) + (f(c1)-0.70);
    end
    totalCost(:,c1) = cumtrapz(1e5*spotSize(:,c1)) + 0.1*(f(c1)-0.70);
end

% normalize time by focal length to account for integration time difference
% considering required SNR
t_norm1 = t/f(1);
t_norm2 = t/f(2);

% runningCost = cumtrapz(cost);
icCross = crossing(spotSize(:,1)-spotSize(:,2));
tcCross = crossing(totalCost(:,1)-totalCost(:,2));

%% Plots
set(0,'DefaultLineLineWidth',2.5)
set(0,'DefaultAxesLineWidth',2.5)
set(0,'DefaultPatchLineWidth',2.5)
set(0,'DefaultTextFontSize',14)
set(0,'DefaultAxesFontSize',14)

fig1 = figure(1);
subplot(211)
plot(t,w,t,u)
ylabel('Torque (N-m)');
% suptitle('Linear Simulation')
legend('Disturbance','Control','Location','NorthEast');
subplot(212)
plot(t,thetaZ)
ylabel('Angle Error (rad)'); xlabel('Time (s)');
% legend('Angle Error','Location','NorthEast');

fig2 = figure(2);
% subplot(311)
% plot(t,u)
% title('Control')
% legend('input','Location','EastOutside')
% ylabel('Torque (N-m)');
subplot(211)
plot(t,spotSize);
hxic = gca;
for cnt = 1:length(icCross)
    x = mean([t(icCross(cnt)) t(icCross(cnt)+1)]);
    line([x x], get(hxic,'YLim'));
end
title('Spot Size')
% xlim([1 6]); ylim([0 3]);
legend('0.70','0.80','Location','SouthEast');
ylabel('Radius (m)')
subplot(212)
plot(t,totalCost);
% plot(t_norm1,runningCost(:,1),t_norm2,runningCost(:,2));
hxrc = gca;
for cnt = 1:length(tcCross)
    x = mean([t(tcCross(cnt)) t(tcCross(cnt)+1)]);
    line([x x], get(hxrc,'YLim'));
end
title('Total Cost')
% xlim([1 6]); ylim([0 5]);
legend('0.70','0.80','Location','SouthEast');
xlabel('Time (s)');

fig3 = figure(3);
plot(t,totalCost);
% plot(t_norm1,runningCost(:,1),t_norm2,runningCost(:,2));
hxrc = gca;
for cnt = 1:length(tcCross)
    x = mean([t(tcCross(cnt)) t(tcCross(cnt)+1)]);
    line([x x], get(hxrc,'YLim'));
end
title('Total Cost')
xlim([0 0.2]); ylim([0 0.5]);
legend('0.70','0.80','Location','West');
xlabel('Time (s)');

% figure(4)
% plot(t,totalCost(:,1)-totalCost(:,2))
% title('Total Cost Difference (0.70 - 0.80)')
% xlabel('Time (s)');



print(fig1,'-depsc','dynamicsPlot')
print(fig2,'-depsc','costPlots')
print(fig3,'-depsc','zoomCostPlot')

%% Conclusions
% In this example, the cost integrand is the "spot size" plus focal length 
% of the mirror. These are related to image quality and field-of-view,
% respectively. The plots show the lowest cost is a function of the focal
% length dependent on the dynamics of the controller. When the disturbance
% is significant, the cost from the spot size is greater. More
% investigation is necessary to deterimine which focal length is optimal.
% The controller parameters might also affect the optimal focal length.
% This is why it is necessary to build an optimal control vector for
% comparison rather than design a feedback system, unless we use an optimal
% feedback controller. In this case, I assume the difference between the
% inertia and size of the different focal length mirrors is negligible,
% though the distance between the mirror and the detector plane changes to
% match the focal length.

% an integration time of 1 second is reasonable
% (http://rsaa.anu.edu.au/files/taros_manual.pdf)