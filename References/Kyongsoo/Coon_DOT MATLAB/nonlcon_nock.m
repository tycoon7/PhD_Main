%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Eigenstructure Assignment Technique
% Planar model control by 3rd and 4th mode change
% Constraint Function File
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq]=nonlcon_nock(x,Ac,Bc,Cc,max_k)
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
for jj=1:2:8
    VWN(:,jj)=real(VWN(:,jj));
    VWN(:,jj+1)=imag(VWN(:,jj+1));
end

% required gain
Kbar=-VWN(9:12,:)*pinv(Cc*VWN(1:8,:));
for jj=1:8
    temp(4*(jj-1)+1:4*(jj-1)+4,1)=Kbar(:,jj);
end

temp=abs(temp);
c=[temp-max_k*ones(32,1);
x];
ceq=[];

