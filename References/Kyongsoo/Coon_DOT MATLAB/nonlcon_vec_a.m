%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Planar model control By Eigenstructure Assignment
% Constraint Function file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq]=nonlcon_vec_a(x,Ac,Bc,Cc,max_k,vd,lamd,Trans)

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

% required gain
Kbar=-VWN(9:12,:)*pinv(Cc*VWN(1:8,:));
m_k=max(max(abs(Kbar)));
c=[m_k-max_k;
-x;
x-ones(2,1)];
ceq=[];