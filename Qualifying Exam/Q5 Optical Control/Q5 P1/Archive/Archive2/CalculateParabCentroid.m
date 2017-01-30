% calculate the centroid of a parabolic arc

clear; close all; clc;

% give information on the conic section
f_val   =  1.0;     % (m)
rho_val =  0.1;     % (kg/m) length density of the mirror
x1_val  = -0.4;     % (m)
x2_val  = -0.2;     % (m)

syms x f rho x1 x2

% parabola equation
y = x^2/(4*f);          % (m)
dydx = x/(2*f);         % (-)
ds = sqrt(1+dydx^2);    % (m) differential arc length

% calculate mass
M_ind = rho*int(ds,x);
M_def = subs(M_ind,x,x2) - subs(M_ind,x,x1);
M = eval(subs(M_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% calculate moment about x-axis
Mx_ind = rho*int(y*ds,x);
Mx_def = subs(Mx_ind,x,x2) - subs(Mx_ind,x,x1);
Mx = eval(subs(Mx_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% calculate moment about y-axis
My_ind = rho*int(x*ds,x);
My_def = subs(My_ind,x,x2) - subs(My_ind,x,x1);
My = eval(subs(My_def,[f rho x1 x2],[f_val rho_val x1_val x2_val]));

% calculate centroid coordinates
x_c = My/M;
y_c = Mx/M;

%% Plot stuff for sanity check

X = linspace(x1_val,x2_val);
Y = subs(y,f,f_val);
Y = eval(subs(Y,x,X));

figure()
hold on
plot(X,Y)
plot(x_c,y_c,'b+','MarkerSize',15)
xlim([-0.5 0.1])
axis equal
hold off
