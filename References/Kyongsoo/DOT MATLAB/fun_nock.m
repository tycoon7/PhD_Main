%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eigenstructure Assignment Technique
% Planar model control by 3rd and 4th mode change
% Objective Function File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=fun_nock(x,Ac,Bc,Cc,max_k)
[vd,lamd]=eig(Ac);lamd=diag(lamd);
wn5=abs(lamd(5)); %Natural Frequency (Mode 4)
wn7=abs(lamd(7)); %Natural Frequency (Mode 3)
lamd(5)=x(1)+sqrt(wn5^2-x(1)^2)*i;
lamd(6)=x(1)-sqrt(wn5^2-x(1)^2)*i;
lamd(7)=x(2)+sqrt(wn7^2-x(2)^2)*i;
lamd(8)=x(2)-sqrt(wn7^2-x(2)^2)*i;

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
E=eig(Ac-Bc*Kbar*Cc);E=real(E);I=find(E>0);s=size(I);
CON=cond(Ac-Bc*Kbar*Cc);
if CON<10^9 & s(1)==0
    Sw=diag([1 1 1 1])*2*10^-7; % Assumed White Gaussian noise Intensity
    %closed loop root mean square(lqr)
    S_x=lyap(Ac-Bc*Kbar*Cc,Bc*Sw*Bc'); %Variance of x matrix
    S_y=Cc(1:4,:)*S_x*Cc(1:4,:)'; %Variance of y matrix
    RMS_cl=diag(sqrt(S_y)); %Close loop Root mean square of x
    P=diag([10^4,10^4,1,1]); %Penalty matrix
    f=RMS_cl'*P*RMS_cl*10^10; %Objective Function
else
    f=10^20;
end