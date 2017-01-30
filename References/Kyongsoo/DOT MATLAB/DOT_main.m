%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eigenstructure Assignment Technique
% Planar model control by changing 3rd and 4th mode
% Main script file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%driving the equations of motion for 2 mirrors and a base
clc;clear;close all;
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
al_0=pi/4-phi/2;

%distance matrix
DM=[1 0 -1 -l*cos(al_0) 0 (L+l*cos(al_0));
1 0 -1 l*cos(al_0) 0 (L- l*cos(al_0));
0 1 -1 0 -l*cos(al_0) -(L- l*cos(al_0));
0 1 -1 0 l*cos(al_0) -(L+l*cos(al_0))];

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
    -sin(phi) 0 sin(phi) tan(phi)*sqrt(L^2+h^2) 0 -(h+hm)*cos(phi) ; %mirror 1 wavelength error
    0 -sin(phi) sin(phi) 0 -tan(phi)*sqrt(L^2+h^2) (h+hm)*cos(phi)]; %mirror 2 wavelength error

C=[Cm zeros(4,6);
   zeros(4,6) Cm];

D=zeros(8,4);

[vo do]=eig(A);
do=diag(do);

%Controllable system Identification
[ABAR,BBAR,CBAR,T,KK] =ctrbf(A,B,C);
Ac=ABAR(5:12,5:12);
Bc=BBAR(5:12,:);
Cc=CBAR(:,5:12);
Gol=ss(Ac,Bc,Cc(1:4,:),D(1:4,:));
[vd,lamd]=eig(Ac);lamd=diag(lamd);
wn5=abs(lamd(5)); %Natural Frequency (Mode 4)
wn7=abs(lamd(7)); %Natural Frequency (Mode 3)
max_k=100;
x0=[real(lamd(5));real(lamd(7))];

%Using Newton's Line Search Method
option=optimset('TolX',1e-12,'TolFun',1e-20,'TolCon',1e-18,'MaxFunEvals',10^6,'MaxIter',10^4,'Display','iter');
[x,fval1,exitflag,output] =...
    fmincon(@fun_nock,x0,[],[],[],[],[],[],@nonlcon_nock,option,Ac,Bc,Cc,max_k);

%Overlapping Desired Eigenvalue
lamd(5)=x(1)+sqrt(wn5^2-x(1)^2)*i;
lamd(6)=x(1)-sqrt(wn5^2-x(1)^2)*i;
lamd(7)=x(2)+sqrt(wn7^2-x(2)^2)*i;
lamd(8)=x(2)-sqrt(wn7^2-x(2)^2)*i;

%Eigenstructure Assignment Technique
for jj=1:8
    N=[lamd(jj)*eye(8)-Ac, -Bc ,zeros(8,8);
    zeros(4,12) ,Bc';
    eye(8) ,zeros(8,4), (lamd(jj)*eye(8)-Ac)'];
    VWN(:,jj)=inv(N)*[zeros(12,1);vd(:,jj)];
end
VWN;

for jj=1:2:8
    VWN(:,jj)=real(VWN(:,jj));
    VWN(:,jj+1)=imag(VWN(:,jj+1));
end
VWN;

% required gain
Kbar=-VWN(9:12,:)*pinv(Cc*VWN(1:8,:));
Sw=diag([1 1 1 1])*2*10^-7; % Assumed White Gaussian noise Intensity

%open loop root mean square(lqr)
S_x=lyap(Ac,Bc*Sw*Bc');%Variance of x matrix
S_y=Cc(1:4,:)*S_x*Cc(1:4,:)'; %Variance of y matrix
RMS_cl=diag(sqrt(S_y)); %Close loop Root mean square of x
p_angle=10^4;
P=diag([p_angle,p_angle,1,1]); %Penalty matrix
f_ol=RMS_cl'*P*RMS_cl*10^10;

%closed loop root mean square(lqr)
S_x=lyap(Ac-Bc*Kbar*Cc,Bc*Sw*Bc');%Variance of x matrix
S_y=Cc(1:4,:)*S_x*Cc(1:4,:)'; %Variance of y matrix
RMS_cl=diag(sqrt(S_y)); %Close loop Root mean square of x
p_angle=10^4;
P=diag([p_angle,p_angle,1,1]); %Penalty matrix
f_cl=RMS_cl'*P*RMS_cl*10^10;

%finding Optimum Desired Eigenvector
for ii=1:2:3
lamd(ii+4)
SUB=orth(inv(A-lamd(ii+4)*eye(12))*B);
rank_sub=rank(SUB);
OL=real(vo(7:12,[1,3,5,7]));
OL_b=lamd(ii+4)*OL;
SUB_2=[OL;OL_b];
Rank_sub2=rank(SUB_2);
SUB_total=[SUB SUB_2];
rank_total=rank(SUB_total);
x0=[1;1;1];
option=optimset('TolX',1e-7,'TolFun',1e-7,'TolCon',1e-9,'MaxFunEvals',10^6,'MaxIter',10^4,'Display','iter');
[x,fval,exitflag,output] = ...
    fmincon(@fun_vec_t,x0,[],[],[],[],[],[],@nonlcon_vec,option,SUB_2,C);
x_f=[x;1];
vd_opt(:,ii)=SUB_2*x_f;
end
for kk=1:2:3
vd_opt(:,kk)=vd_opt(:,kk)/norm(vd_opt(:,kk)); %Optimum Eigenvector
end
Trans=T*vd_opt;

%Assigning Proper Ratio between Open-Loop Eigenvector and Optimum Eigenvector
x0=[0.5;0.5];
option=optimset('TolX',1e-7,'TolFun',1e-7,'TolCon',1e-9,'MaxFunEvals',10^6,'MaxIter',10^4,'Display','iter');
[x,fval,exitflag,output] = ...
    fmincon(@fun_vec_a,x0,[],[],[],[],[],[],@nonlcon_vec_a,option,Ac,Bc,Cc,max_k,vd,lamd,Trans);
vd(:,5)=x(1)*vd(:,5)+(1- x(1))*Trans(5:12,1);
vd(:,6)=conj(vd(:,5));
vd(:,7)=x(2)*vd(:,7)+(1- x(2))*Trans(5:12,3);
vd(:,8)=conj(vd(:,7));

for jj=1:8
    N=[lamd(jj)*eye(8)-Ac, -Bc ,zeros(8,8);
    zeros(4,12) ,Bc';
    eye(8) ,zeros(8,4), (lamd(jj)*eye(8)-Ac)'];
    VWN(:,jj)=inv(N)*[zeros(12,1);vd(:,jj)];
end
for jj=1:2:8
    VWN(:,jj)=real(VWN(:,jj));
    VWN(:,jj+1)=imag(VWN(:,jj+1));
end

% required gain matrix
Kbar=-VWN(9:12,:)*pinv(Cc*VWN(1:8,:));
[vc dc]=eig(A-B*Kbar*C);dc=diag(dc);

%closed loop root mean square(lqr)
S_x=lyap(Ac-Bc*Kbar*Cc,Bc*Sw*Bc');%Variance of x matrix
S_y=Cc(1:4,:)*S_x*Cc(1:4,:)'; %Variance of y matrix
RMS_cl=diag(sqrt(S_y)); %Close loop Root mean square of x
f_cl=RMS_cl'*P*RMS_cl*10^10;