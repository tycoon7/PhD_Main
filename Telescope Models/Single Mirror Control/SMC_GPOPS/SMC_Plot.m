%% -----------------------------------------------------------------------%
%------------------------------- Plot Solution ---------------------------%
%-------------------------------------------------------------------------%
figure('position',[100 50 1100 800])
subplot(221)
% plot the states
p1 = plot(solution.phase(1).time, solution.phase(1).state(:,1:2), '-o');
xl = xlabel('time');
yl = ylabel('state');
set(p1,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
h1 = legend('$e$','$\dot{e}$','Location','eastoutside');   
set(h1,'interpreter','latex')
grid on

subplot(222)
% plot the control
p2 = plot(solution.phase(1).time,solution.phase(1).control,'-o');
xl = xlabel('time');
yl = ylabel('control');
set(p2,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
h2 = legend('$u$','Location','eastoutside');
set(h2,'interpreter','latex')
grid on

subplot(223)
% plot the Control Torque
% pp = plot(solution.phase(1).time,sin(10*solution.phase(1).time),'-o');
w = 1*sin(10*linspace(t0,tf,1000)+auxdata.ps);
p3 = plot(linspace(t0,tf,1000),w,solution.phase.time,solution.phase.state(:,3));
xl = xlabel('time');
yl = ylabel('Torque');
set(p3,'LineWidth',1.25,'MarkerSize',8);
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontSize',16);
h3 = legend('$w$','$\tau_u$','Location','eastoutside');
set(h3,'interpreter','latex')
grid on

% calculate the spot size wrt time using the final states
e = solution.phase.state(:,1);
max_e = max(e);
spotSize = zeros(length(e),1);
% for i = 1:length(e)
%     spotSize(i,1) = SingleMirror_SpotSize(e(i),auxdata.f);
% end
% max_ss = max(spotSize);


hp4 = subplot(224);
focalLength = num2str(auxdata.f,'%4.2f');
fString = strcat('Focal Length = ',focalLength);
text(0.25,0.8,fString);
inttext='Integral = ';
int = num2str(trapz(abs(e)),4);
intstring = strcat(inttext,int);
text(0.25,0.6,intstring);
% maxSStext = 'Maximum Spot Size = ';
% maxSS = num2str(max_ss);
% SSstring = strcat(maxSStext,maxSS);
% text(0.25,0.4,SSstring);
maxEtext = 'Maximum Angular Error = ';
maxE = num2str(max_e);
Estring = strcat(maxEtext,maxE);
text(0.25,0.4,Estring);
set(hp4,'visible','off');

% % plot the cost integrand wrt time
% p4 = plot(solution.phase.time,spotSize);
% xl = xlabel('time');
% yl = ylabel('RMS Spot Size');
% % yl = ylabel('Angle Error');
% set(p4,'LineWidth',1.25,'MarkerSize',8);
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontSize',16);
% % h4 = legend('$\tau$','Location','eastoutside');
% % set(h4,'interpreter','latex')
% grid on

filename = strcat('f_',focalLength,'.fig');
savefig(filename);
