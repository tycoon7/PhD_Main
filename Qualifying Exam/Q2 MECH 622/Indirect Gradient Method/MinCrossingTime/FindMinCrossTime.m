% find the velocity of the river to give shortest crossing time
clear; close all; clc;

Vr_knots = 0:0.1:2;          % (knots) river velocities vector
Vr_ms = Vr_knots*0.514;     % (m/s) river velocities
u = zeros(51,1);
save('u0_guess.mat','u');

for i = 1:length(Vr_knots)
    I = strcat('Vr_knots = ', num2str(Vr_knots((i))));
    disp(I)
    MECH622_Qual;
    crossingTime(i) = tf;
%     title = strcat('CrossingPath_',num2str(Vr_knots(i)),'_knots');
%     savefig(title);
%     save(title);
    save('u0_guess.mat','u');
end

[minT, index] = min(crossingTime);

%% Plots
figure()
hold on
plot(Vr_knots,crossingTime)
plot(Vr_knots(index),minT,'ro')
hold off
xlabel('Max River Velocity (knots)');
ylabel('Min Time to Cross (s)');
title({'Minimum Cross Time is 120.9 sec'; ...
       'when Vmax = 0.7 knots'});

% Title1 = strcat({'Minimum Cross Time is '}, num2str(minT), {' sec'});
% Title2 = strcat({'when Vmax = '},num2str(Vr_knots(index)),{' knots'});
% title({ Title1; Title2});

