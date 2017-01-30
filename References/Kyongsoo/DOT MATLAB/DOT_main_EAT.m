%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AFRL Model control By Eigenstructure Assignment Technique
% Main script file
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clc;clear;close all;
% Model Reduction
load tuned_best_fit_0L
sys=d2c(fit_sys);
n=size(sys(:,1));
h=minreal(sys);
[hb,g]=balreal(h);
index=21:79;
hdel=modred(hb,index,'del');
figure(1);sigma(hdel,{200,2000});grid on
[A,B,C,D]=ssdata(hdel);
%Setting Open Loop Value as desired one.
[vd lamd]=eig(A);lamd=diag(lamd);
[xx0,I]=sort(abs(lamd));lamdd=lamd;vdd=vd;
%Arranging Eigenvalue by increasing order
for kk=1:20
lamd(kk)=lamdd(I(kk));
vd(:,kk)=vdd(:,I(kk));
end
for kk=1:9
r(kk,1)=real(lamd(2*kk));
Wn(kk,1)=abs(lamd(2*kk));
end
max_k=1 %Max allowable gain matrix element.
%Open loop root mean square(lqr)
load XX_modred % Contains variable for the White Gaussian noise
intensity (x)
Sw=diag(x); % White Gaussian noise intensity
S_x=lyap(A,B*Sw*B'); %Variance of x matrix
S_y=C*S_x*C'; %Variance of y matrix
RMS1=diag(sqrt(S_y)) %Closed loop Root mean square of x
f=RMS1'*RMS1/10^4 %Objective Function
%Finding optimum eigenvalue assignment with open loop eigenvector for desired one.
x0=[r]
option=optimset('TolX',1e-9,'TolF un',1e-9,'TolCon',1e-9,'MaxFunEvals',10^6,'MaxIter',40^3,'Display','iter');
[x,fval1,exitflag,output] =
fmincon(@fun_e_val,x0,[],[],[],[],[],[],@nonlcon_e_val,option,A,B,C,max_k,lamd,vd,
Wn);
save result_e_val x
for kk=1:9
lamd(2*kk-1)=x(kk,1)+sqrt(Wn(kk)^2-x(kk,1)^2)*i;
lamd(2*kk)=conj(lamd(2*kk-1));
end
% Checking the result of assigning the eigenvalue.
P=eye(20); %eigenvector penalty matrix
for jj=1:20
N=[lamd(jj)*eye(20)-A, -B ,zeros(20);
zeros(9,20+9) ,B';
P ,zeros(20,9), (lamd(jj)*eye(20)-A)'];
VWN(:,jj)=inv(N)*[zeros(29,1);P*vd(:,jj)];
end
for jj=1:2:20
VWN(:,jj)=real(VWN(:,jj));
VWN(:,jj+1)=imag(VWN(:,jj+1));
end
Kbar=-VWN(20+1:20+9,:)*pinv(eye(20)*VWN(1:20,:)); % required gain matrix with
full state-feedback
G_cl=ss(A-B*Kbar,B,C,D);
figure(2);sigma(G_cl,{200,2000});grid on
%Closed loop root mean square(lqr)
load XX_modred % Contains variable for the White Gaussian noise
Sw=diag(x); % White Gaussian noise intensity
S_x=lyap(A-B*Kbar,B*Sw*B'); %Variance of x matrix
S_y=C*S_x*C'; %Variance of y matrix
RMS=diag(sqrt(S_y)) %Closed loop Root mean square of x
f=RMS'*RMS/10^11
MAX_K=max(max(abs(Kbar))) %Objective Function
save E_VAL_result