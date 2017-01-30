%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Planar model control By L Q R Method
% Constraint Function file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq]=nonlcon_lqr(x,Ac,Bc,Cc,max_k)
Q=diag([x])*10^1; % State Penalty
R=diag([1,1,1,1]*.1);% Control Penlaty
K=lqr(Ac,Bc,Q,R);
Kbar=K*Cc^-1;
for jj=1:8
temp(4*(jj-1)+1:4*(jj-1)+4,1)=Kbar(:,jj);
end
temp=abs(temp);
c=[temp-max_k*ones(32,1);
-x+0.1*ones(8,1)];
ceq=[];