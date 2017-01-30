%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Planar model control By L Q R Method
% Main script file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%driving the equations of motion for 2 mirrors and a base
clc
clear
close all
I1= 1.5;
I2= 1.5;
Ib= 5;
m1= 1;
m2= 1;
mb= 50;
h= 4;
hm= 0.1;
L= 1;
l= 0.1;
k= 2.5*10^6;
c= 100;
%Initial Set up
phi=atan(h/L);
al_0=pi/4-phi/2 ;
DM=[1 0 -1 -l*cos(al_0) 0 (L+l*cos(al_0));
1 0 -1 l*cos(al_0) 0 (L- l*cos(al_0));
0 1 -1 0 -l*cos(al_0) -(L- l*cos(al_0));
0 1 -1 0 l*cos(al_0) -(L+l*cos(al_0))];%distance matrix
DDM=[DM zeros(4,6);
zeros(4,6) DM];
ACM=[k 0 0 0 c 0 0 0;
0 k 0 0 0 c 0 0;
0 0 k 0 0 0 c 0;
0 0 0 k 0 0 0 c];

AB=[1/m1*[-1 -1 0 0];
1/m2*[0 0 -1 -1];
1/mb*[1 1 1 1];
1/I1*l*cos(al_0)*[1 -1 0 0];
1/I2*l*cos(al_0)*[0 0 1 -1];
1/Ib*[(L+l*cos(al_0))*[-1 0 0 1]+(L-l*cos(al_0))*[0 -1 1 0]]];
A=[zeros(6) eye(6)
AB*ACM*DDM];
B=[zeros(6,4);
AB];
Cm=[[cos(phi)/sqrt(L^2+h^2), 0, -cos(phi)/sqrt(L^2+h^2), 1, 0, -sin(phi)*(h+hm)/sqrt(L^2+h^2) ; %mirror 1 tilt angle
0, -cos(phi)/sqrt(L^2+h^2), cos(phi)/sqrt(L^2+h^2), 0, 1, -sin(phi)*(h+hm)/sqrt(L^2+h^2) ]; %mirror 2 tilt angle
[-sin(phi) 0 sin(phi) tan(phi)*sqrt(L^2+h^2) 0 -(h+hm)*cos(phi) ; %mirror 1 wavelength error
0 -sin(phi) sin(phi) 0 -tan(phi)*sqrt(L^2+h^2) (h+hm)*cos(phi)]]; %mirror 2 wavelength error
C=[Cm zeros(4,6);
zeros(4,6) Cm];
D=zeros(8,4);
[vo do]=eig(A);do=diag(do);
%Controllable system Identification
[ABAR,BBAR,CBAR,T,KK] =ctrbf(A,B,C);
Ac=ABAR(5:12,5:12);
Bc=BBAR(5:12,:);
Cc=CBAR(:,5:12);
max_k=100;
x=[100;100;100;100;1;.1;.01;.01];
option=optimset('TolX',1e-10,'TolFun',1e-19,'TolCon',1e-10,'MaxFunEvals',10^6,'MaxIter',10^5,'Display','iter');

for kkk=1:1
    load XXX0
    x0=(x*1)/1;
    [x,fval1,exitflag,output] = fmincon(@fun_lqr,x0,[],[],[],[],[],[],@nonlcon_lqr,option,Ac,Bc,Cc,max_k);
    %save XXX0 x;
end

Q=diag([x])*10^1;
R=diag([1,1,1,1]*.1);
K=lqr(Ac,Bc,Q,R);
Kbar=K*Cc^-1;
Sw=diag([1 1 1 1])*2*10^-7; % Assumed White Gaussian noise Intensity
%closed loop root mean square(lqr)
S_x=lyap(Ac-Bc*Kbar*Cc,Bc*Sw*Bc'); %Variance of x matrix
S_y=Cc(1:4,:)*S_x*Cc(1:4,:)'; %Variance of y matrix
RMS_cl=diag(sqrt(S_y)); %Close loop Root mean square of x
P=diag([10^4,10^4,1,1]); %Penalty matrix
f=RMS_cl'*P*RMS_cl*10^10;
