%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% AFRL Model control By Eigenstructure Assignment Technique
% Constraint Function file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [c,ceq]=nonlcon_e_val(x,A,B,C,max_k,lamd,vd,Wn)
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
Kbar=-VWN(20+1:20+9,:)*pinv(eye(20)*VWN(1:20,:)); % required gain matrix
K_max=max(max(abs(Kbar)));
c=[x(:,1);K_max-max_k];
ceq=[];