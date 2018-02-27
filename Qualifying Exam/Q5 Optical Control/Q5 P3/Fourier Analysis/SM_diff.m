function [ diffData ] = SM_diff( run,isCentered,isIdeal )
%SM_DIFF sets up a single mirror and calculates diffraction pattern data
%   Set up a single parabolic mirror to be offset, then add a surface
%   roughness term as a deviation to the gaussian reference sphere (i.e.
%   create the generalized pupil function

fn = 10;                % (-) f-number
D = 25e-3;              % (m) lens diameter
f = fn*D;               % (m) focal length of the lens
L1 = 250e-3;            % (m) source field width
M = 500;                % given number of sample points
dx1=L1/M;
x1=-L1/2:dx1:L1/2-dx1;  % src coords
y1=x1;
z = f;                  % (m) detector plane distance
lambda = 0.5e-6;        % (m) wavelength
k = 2*pi/lambda;        % (1/m) wave number

[X1,Y1]=meshgrid(x1,y1);

% run different scenarios.
% Option #1: centered or offset pupil
% Option #2: deviated or ideal surface
if isCentered
    P = sqrt(X1.^2+Y1.^2) <= D/2;       % pupil function (centered)
else
P = sqrt((X1+D/2).^2+Y1.^2) <= D/2;     % pupil function (shifted in x)
end

if isIdeal
    u1 = P;
else
    W = calcD(X1,Y1,f);          % effective path-length deviation
    gP = P.*exp(1i*k*W);         % generalized pupil function
    u1 = gP;
end

Imin = 151;
Imax = 350;
Ispan = Imin:Imax;
Imid = (Imax-Imin+1)/2 + 1;

% Input field intesity distribution
I1=abs(u1.^2);                      % src irradiance
% limit plot extents
x1_red = x1(Ispan);
y1_red = y1(Ispan);
I1_red = I1(Ispan,Ispan);

% Fraunhofer approximation
[u2,L2]=propFF(u1,L1,lambda,z);
dx2=L2/M;
x2=-L2/2:dx2:L2/2-dx2; %obs ords
y2=x2;
I2=abs(u2.^2);
% limit plot extents
x2_red = x2(Ispan);
y2_red = y2(Ispan);
I2_red = I2(Ispan,Ispan);

%% Plots
if run == 1
    figure(1)
    suptitle({'Focal Plane Diffraction Pattern';''})
end

% title
ax = subplot(3,6,1+6*(run-1));
set(ax,'visible','off')
text(-1,0.5,titles(run),'fontsize',14)
% source
subplot(3,6,2+6*(run-1))
imagesc(x1_red,y1_red,I1_red);
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
if run==1; title({'Input Field';''}); end;
% focal plane intensity diffraction pattern
subplot(3,6,3+6*(run-1))
imagesc(x2_red,y2_red,nthroot(I2_red,5));   %stretch image contrast
axis square; axis xy;
colormap('gray'); xlabel('x (m)'); ylabel('y (m)');
if run==1; title({'Intensity Dist';''}); end;
% cross section of intensity
subplot(3,6,4+6*(run-1))
plot(x2_red,I2(Imid,Ispan));
axis square; xlabel('x (m)'); ylabel('Irradiance');
if run==1; title({'Intensity Cross';''}); end;
% cross section of magnitude
subplot(3,6,5+6*(run-1))
plot(x2_red,abs(u2(Imid,Ispan)));
axis square; xlabel('x (m)'); ylabel('Magnitude');
if run==1; title({'Mag Cross';''}); end;
% cross section of phase
subplot(3,6,6+6*(run-1))
plot(x2_red,unwrap(angle(u2(Imid,Ispan))));
axis square; xlabel('x (m)'); ylabel('Phase (rad)');
if run==1; title({'Phase Cross';''}); end;

end

