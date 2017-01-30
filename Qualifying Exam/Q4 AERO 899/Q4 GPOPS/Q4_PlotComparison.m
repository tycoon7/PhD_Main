%---------------------------------------------------%
% AFIT Qualifying Exam Question #4, 01 Dec 2014     %
% Timothy Coon                                      %
% Advisor: Dr Cobb                                  %
%---------------------------------------------------%

% Q4_PlotComparison.m

x = solution.phase.time;
y = solution.phase.state;
t_cost = solution.phase.integral;
M = length(N_mesh);
C = length(N_cols);
num = C*(i-1) + j;

% figure()
subplot(M,C,num)
plot(x,y,'k-o','linewidth',2)
hold on
Q4_GradientPlot
hold off
title(['Total Cost = ',num2str(t_cost)])
xlim([x0 xf]); ylim([ymin ymax]);
axis equal

if j == 1
    h1 = text(-0.75,0.8,['Num Mesh = ', num2str(i)]);
    set(h1,'rotation',90,'fontsize',14)
end

if i == 1
    h2 = text(0.9,3.65,['Col Pts = ', num2str(j)]);
    set(h2,'fontsize',14)
end
    
if num == M*C
    suptitle('Minimum Energy Path to the Road (No Control Cost)')
end