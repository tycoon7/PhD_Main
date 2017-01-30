% test_calcVvector
% clear all;
close all; clc;

N = 10;
A = pi/2*ones(1,N);
% A = linspace(0,pi/2,N);
% th = -pi/6*ones(1,N);
th = linspace(-pi/4,pi/4,N);
for i = 1:N
    V(:,i) = calcVvector(A(i),th(i));
end

figure(1)
plot(V(1,:),'linewidth',2)
hold on
plot(V(2,:),'linewidth',2)
hold off
legend('Vx','Vy')