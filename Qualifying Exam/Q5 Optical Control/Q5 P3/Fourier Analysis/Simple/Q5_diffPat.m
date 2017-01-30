% Lens example, page 97
% also, start of problem 6.4
clear; close all; clc;

fn = 10;                % (-) f-number
D = 25e-3;              % (m) lens diameter
f = fn*D;               % (m) focal length of the lens
L1 = 250e-3;            % (m) source field width
M = 250;                % given number of sample points
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1;  % src coords
y1=x1;
z = f;                  % (m) detector plane distance
lambda = 0.5e-6;        % (m) wavelength
k = 2*pi/lambda;        % (1/m) wave number

[X1,Y1]=meshgrid(x1,y1);
% P = sqrt(X1.^2+Y1.^2) <= D/2;           % pupil function (centered)
P = sqrt((X1+D/2).^2+Y1.^2) <= D/2;     % pupil function (shifted in x)
W = calcD(X1,Y1,f);                            % effective path-length deviation
gP = P.*exp(1i*k*W);                     % generalized pupil function
u1 = gP;

I1=abs(u1.^2);                      % src irradiance

figure(1)
imagesc(x1,y1,I1);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
title('z = 0 m');

[u2,L2]=propFF(u1,L1,lambda,z);

dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; %obs ords
y2=x2;

I2=abs(u2.^2);

figure(2)
suptitle('Focal Plane Diffraction Pattern')
subplot(221)
imagesc(x2,y2,nthroot(I2,4));   %stretch image contrast
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
% title(['z= ',num2str(z),' m']);
%
subplot(222) %irradiance profile
plot(x2,I2(M/2+1,:));
axis square; xlabel('x (m)'); ylabel('Irradiance');
% title(['z= ',num2str(z),' m']);
%
subplot(223) %plot obs field mag
plot(x2,abs(u2(M/2+1,:)));
axis square; xlabel('x (m)'); ylabel('Magnitude');
% title(['z= ',num2str(z),' m']);
%
subplot(224) %plot obs field phase
plot(x2,unwrap(angle(u2(M/2+1,:))));
axis square; xlabel('x (m)'); ylabel('Phase (rad)');
% title(['z= ',num2str(z),' m']);