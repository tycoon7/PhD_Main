function [Pc,P1,P2,M,I] = calcConicCentroid(x1_val,x2_val,f_val,rho_val)
% calculate the centroid of a 2d parabolic arc
% can be easily modified for any conic by changing the equation
% look into making this more dynamic by including the eccentricity
% Pc = vector from vertex to CoM
% P1 = vector from vertex to point 1
% P2 = vector from vertex to point 2
% M  = mass of mirror
% I  = moment of inertia of mirror

% give information on the conic section
% x1_val  = -0.4;     % (m)
% x2_val  = -0.2;     % (m)
% f_val   =  1.0;     % (m)
% rho_val =  0.1;     % (kg/m) length density of the mirror

syms x f rho x1 x2

% parabola equation
y = x^2/(4*f);          % (m)
dydx = x/(2*f);         % (-)
ds = sqrt(1+dydx^2);    % (m) differential arc length

% calculate mass
M_ind = rho*int(ds,x);
M_def = subs(M_ind,x,x2) - subs(M_ind,x,x1);
M = eval(subs(M_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% calculate first moment about x-axis
Mx_ind = rho*int(y*ds,x);
Mx_def = subs(Mx_ind,x,x2) - subs(Mx_ind,x,x1);
Mx = eval(subs(Mx_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% calculate first moment about y-axis
My_ind = rho*int(x*ds,x);
My_def = subs(My_ind,x,x2) - subs(My_ind,x,x1);
My = eval(subs(My_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% calculate centroid coordinates
x_c = My/M;
y_c = Mx/M;

% calculate the mass moment of inertia about the centroid
rx = x - x_c;
ry = y - y_c;
r = sqrt(rx^2 + ry^2);
I_ind = rho*int(r^2,x);
I_def = subs(I_ind,x,x2) - subs(I_ind,x,x1);
I = eval(subs(I_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% ouput values
Pc = [x_c; y_c];
P1 = [x1_val;eval(subs(y,[f,x],[f_val,x1_val]))];
P2 = [x2_val;eval(subs(y,[f,x],[f_val,x2_val]))];

%% Plot stuff for sanity check

% X = linspace(x1_val,x2_val);
% Y = subs(y,f,f_val);
% Y = eval(subs(Y,x,X));
% Pc = [x_c; y_c];
% P1 = [X(1);Y(1)];
% P2 = [X(end);Y(end)];
% 
% figure()
% hold on
% plot(X,Y)
% plot(x_c,y_c,'b+','MarkerSize',15)
% xlim([-0.5 0.1])
% axis equal
% hold off
