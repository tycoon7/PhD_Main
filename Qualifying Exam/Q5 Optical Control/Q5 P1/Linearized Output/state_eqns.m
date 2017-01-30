function zdot = state_eqns(t,z,Dt,D)
% state equations for qualifying exam question 5, problem 1

global I m c k rho1 rho2 phi1 phi2

zdot = zeros(6,1);

% Control Input
% Ux = U(:,1);
% Ux = interp1(Ut,Ux,t);
% Uy = U(:,2);
% Uy = interp1(Ut,Uy,t);

% Disturbance Input
Dx = D(:,1);
Dx = interp1(Dt,Dx,t);
Dy = D(:,2);
Dy = interp1(Dt,Dy,t);

% first-Order State Space Model
x      = z(1);
xd     = z(2);
y      = z(3);
yd     = z(4);
alpha  = z(5);
alphad = z(6);

x1  = x  + rho1(2)*alpha*cos(phi1);
x1d = xd + rho1(2)*alphad*cos(phi1);
x2  = x  + rho2(2)*alpha*cos(phi2);
x2d = xd + rho2(2)*alphad*cos(phi2);
y1  = y  + rho1(1)*alpha*sin(phi1);
y1d = yd + rho1(1)*alphad*sin(phi1);
y2  = y  + rho2(1)*alpha*sin(phi2);
y2d = yd + rho2(1)*alphad*sin(phi2);

Sx1 = c*x1d + k*x1;
Sx2 = c*x2d + k*x2;
Sy1 = c*y1d + k*y1;
Sy2 = c*y2d - k*y2;

zdot(1) = xd;
zdot(2) = (1/m)*(-Sx1 - Sx2 + Dx);
zdot(3) = yd;
zdot(4) = (1/m)*(-Sy1 - Sy2 + Dy);
zdot(5) = alphad;
zdot(6) = (1/I)*(rho2(2)*(-Sy2 + Dx)...
               + rho2(1)*(-Sy2 + Dy)...
               + rho1(2)*Sx1 + rho1(1)*Sy1);
end