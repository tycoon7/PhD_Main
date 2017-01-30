% Realize the notion of the surface roughness function on a parabolic arc
clear; close all; clc;

syms lsym
f = 0.2;      % focal length of parabolic mirror

%% ideal parabolic arc: x^2 = 4*f*y
x0 = 0; xf = -0.25;
% t0 = x0/(2*f); tf = xf/(2*f);
t0 = 0; tf = (xf-x0)/(2*f);
t = linspace(t0,tf,100);
% xi = 2*f*(t + (x0/(2*f)));                     % ideal surface x-pos
xi = 2*f*t + x0;
dxi = 2*f*ones(1,length(t)); 
yi = f*t.^2;                    % ideal surface y-pos
dyi = 2*f*t;
Pi = [xi; yi];
dPi = [dxi; dyi];

%% calculate total arclength
l2 = f*(asinh(tf) + tf*sqrt(tf^2+1));   % to end point
l1 = f*(asinh(t0) + t0*sqrt(t0^2+1));   % to start point
l0 = l1 - l2;

%% calculate the deviation from the ideal surface
l = f*(asinh(t) + t.*sqrt(t.^2+1));
d = 0;              % surface roughness term
dd = 0;             % derivative of deviation
di = [1; 1; 1; 1; 1]*1.0e-5;  % weights of deviation (~ wavelengths)
% build the symbolic function
for i = 1:length(di)
    d = d + di(i)*sin(pi*i*lsym/l0);
    dd = dd + di(i)*cos(pi*i*lsym/l0)*pi*i/l0;
end
% evaluate the function
for j = 1:length(l)
    D(j,1) = subs(d,lsym,l(j));
    DD(j,1) = subs(dd,lsym,l(j));
end
D = eval(D);
DD = eval(DD);

%% calculate ideal surface normals
% This negative reciprocal only works on left side of y-axis 
normdPi = sqrt(sum(abs(dPi).^2,1));       % two-norm of each column
ni = [-dyi./normdPi; dxi./normdPi]; 

%% calculate actual surface coords
x = xi + ni(1,:).*D.';
dx = dxi + ni(1,:).*DD.';
y = yi + ni(2,:).*D.';
dy = dyi + ni(2,:).*DD.';
P = [x; y];
dP = [dx; dy];

%% calculate actual surface normals
normdP = sqrt(sum(abs(dP).^2,1));       % two-norm of each column
n = [-dy./normdP; dx./normdP]; 

%% find focal plane intersection points
% create index vector to pick out intersection points, generate reflected 
% ray direction vectors, then line functions, then find intersection

i0_hat = [0; -1];               % direction of incident rays
index = [1:10:100 100];         % index, include start and endpoint
startPnt = [x(index); 0.25*ones(1,length(index))];
R = zeros(2,2,length(index));
r_hat = zeros(2,length(index));
for i = 1:length(index)
    j = index(i);
    R(:,:,i) = eye(2) - 2*n(:,j)*n(:,j).';  % reflection dyadic
    r_hat(:,i) = R(:,:,i)*i0_hat;
end
m = r_hat(2,:)./r_hat(1,:);
m(m==Inf) = 0;
b = y(index)-m.*x(index);
syms x1
y_ray = m.*x1+b;
for i = 1:length(index)
    if y_ray(i) == 0
        x2(i) = 0;
    else
        x2(i) = solve(y_ray(i)-f,x1); % solve is fussy, zeros are doubles
    end
end
x2 = eval(x2);      % messes up with zeros, sometimes

%% calculate RMS spot size radius
RMS_radius = max(x2) - min(x2);

%% plots
scrsz = get(0,'ScreenSize');
scrwidth = scrsz(3);
scrheight = scrsz(4);
pos1 = [0 0 scrwidth/2 scrheight];
pos2 = [scrwidth/2 0 scrwidth/2 scrheight];

figure('Position',pos1)
subplot(211)
plot(xi,yi,x,y,'--r')
hold on
arrow(P(:,index).',P(:,index).'+0.04*n(:,index).')
hold off
axis equal
ylabel('y'); xlabel('x');
title('Ideal and Actual Surfaces (With Normal Vectors)');
subplot(212)
plot(x,D)
ylabel('Deviation Along Normal'); xlabel('x-Position');
title(strcat('Surface Roughness. Weights = [1; 1; 1; 1; 1]*1.0e-5'));

figure('Position',pos2)
plot(xi,yi,x,y,'--r')
hold on
plot(x(index),y(index),'*',x2,f,'b+')
arrow([x(index)',y(index)'],[x2' f*ones(length(index),1)])
arrow(startPnt,[x(index)',y(index)'])
hold off
axis equal; axis square;
ylabel('y'); xlabel('x');
title(strcat('Rough Parabolic Mirror. RMS Spot Size = ',...
                                                     num2str(RMS_radius)));
