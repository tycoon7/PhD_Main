%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AFRL Model control By Eigenstructure Assignment Technique
% Objective Function file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=fun_e_val(x,A,B,C,max_k,lamd,vd,Wn)

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
Kbar=-VWN(20+1:20+9,:)*pinv(eye(20)*VWN(1:20,:)); % required gain
E=eig(A-B*Kbar);E=real(E);I=find(E>0);s=size(I);
CON=cond(A-B*Kbar);
if CON<10^9 & s(1)==0
    load XX_modred % Contains variable for the White Gaussian noise intensity
    Sw=diag(x); % White Gaussian noise intensity
    %Closed loop root mean square(lqr)
    S_x=lyap(A-B*Kbar,B*Sw*B'); %Variance of x matrix
    S_y=C*S_x*C'; %Variance of y matrix
    RMS=diag(sqrt(S_y)); %Closed loop Root mean square of x
    f=RMS'*RMS/10^4;
else
    f=10^20;
end