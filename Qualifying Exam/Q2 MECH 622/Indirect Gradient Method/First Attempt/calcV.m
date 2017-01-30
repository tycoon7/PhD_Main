function F = calcV(x,u,Vriver)
% calcV holds the set of equations to calculate relative velocity of the
% boat and other stuff using fsolve

global Vwx Vwy

th = u(1);
A = u(2);
psi = deg2rad(120)-A;
W = sqrt((Vwx-Vriver)^2 + Vwy^2);

% guess a value of mu (will depend on boat properties)
mu = 1;         % coefficients

V = x(1);
Wr = x(2);
alpha = x(3);

F = [mu^2*Wr^2*sin(alpha)*sin(th) - V^2;
     W*sin(psi) - Wr*sin(alpha+th);
     V^2 + W^2 - 2*V*W*cos(psi) - Wr^2];

