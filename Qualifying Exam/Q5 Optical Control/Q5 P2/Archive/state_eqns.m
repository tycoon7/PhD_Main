function zdot = state_eqns(t,z,Dt,D,Ut,U)
% state equations for qualifying exam question 5, problem 2

global R I m c k

zdot = zeros(6,1);

Ux = U(:,1);
Ux = interp1(Ut,Ux,t);
Uy = U(:,2);
Uy = interp1(Ut,Uy,t);

Dx = D(:,1);
Dx = interp1(Dt,Dx,t);
Dy = D(:,2);
Dy = interp1(Dt,Dy,t);

x      = z(1);
xd     = z(2);
y      = z(3);
yd     = z(4);
alpha  = z(5);
alphad = z(6);

rho1 = [-(R-(4*R/(3*pi))); 4*R/(3*pi)];
rho2 = [4*R/(3*pi); -(R-(4*R/(3*pi)))];
phi1 = atan((R-(4*R/(3*pi)))/(4*R/(3*pi)));
phi2 = atan((4*R/(3*pi))/(R-(4*R/(3*pi))));

x1  = x  + rho1(2)*alpha*cos(phi1);
x1d = xd + rho1(2)*alphad*cos(phi1);
x2  = x  + rho2(2)*alpha*cos(phi2);
x2d = xd + rho2(2)*alphad*cos(phi2);
y1  = y  + rho1(1)*alpha*sin(phi1);
y1d = yd + rho1(1)*alphad*sin(phi1);
y2  = y  + rho2(1)*alpha*sin(phi2);
y2d = yd + rho2(1)*alphad*sin(phi2);
x0  = x1 + Ux;
x0d = x1d;
y0  = y1 + Uy ;
y0d = y1d;

zdot(1) = xd;
zdot(2) = (1/m)*(-c*x0d - k*x0 - c*x2d - k*x2 + Dx);
zdot(3) = yd;
zdot(4) = (1/m)*(-c*y0d - k*y0 - c*y2d - k*y2 + Dy);
zdot(5) = alphad;
zdot(6) = (1/I)*((4*R/(3*pi))*(-c*x0d - k*x0 - c*y2d - k*y2 + Dy)...
               + (R-(4*R/(3*pi)))*(-c*x2d - k*x2 + Dx - c*y0d - k*y0));

% zdot(1) = xd;
% zdot(2) = (1/m)*(-c*x1d - k*x1 - c*x2d - k*x2 + Dx);
% zdot(3) = yd;
% zdot(4) = (1/m)*(-c*y1d - k*y1 - c*y2d - k*y2 + Dy);
% zdot(5) = alphad;
% zdot(6) = (1/I)*((4*R/(3*pi))*(-c*x1d - k*x1 - c*y2d - k*y2 + Dy)...
%                + (R-(4*R/(3*pi)))*(-c*x2d - k*x2 + Dx - c*y1d - k*y1));
           
           
