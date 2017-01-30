%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Planar model control By L Q R Method
% Objective Function file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f=fun_lqr(x,Ac,Bc,Cc,max_k)
Q=diag([x])*10^1; % State Penalty
R=diag([1,1,1,1]*.1);% Control Penlaty
K=lqr(Ac,Bc,Q,R);
Kbar=K*Cc^-1;
Sw=diag([1 1 1 1])*2*10^-7; % Assumed White Gaussian noise Intensity
%closed loop root mean square(lqr)
S_x=lyap(Ac-Bc*Kbar*Cc,Bc*Sw*Bc');%Variance of x matrix
S_y=Cc(1:4,:)*S_x*Cc(1:4,:)'; %Variance of y matrix
RMS_cl=diag(sqrt(S_y)); %Close loop Root mean square of x
P=diag([10^4,10^4,1,1]); %Penalty matrix
f=RMS_cl'*P*RMS_cl*10^10;