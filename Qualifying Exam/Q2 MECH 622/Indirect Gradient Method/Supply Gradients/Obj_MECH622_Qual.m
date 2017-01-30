function [phi,Huphid,s]=Obj_MECH622_Qual(ut,s0,name,dims)
% ut = vector of controls
% s0 = initial state vector
% name = name of function to propagate and other calcs
% dims = vector of pertinent dimensions

N=dims(1);
ns=dims(2);
nc=dims(3);
ncons=dims(4);
nt1=ncons+1;
N1=N+1;
s=zeros(ns,N1);
la=zeros(ns,nt1);
Hu=zeros(nt1,nc,N);
s(:,1)=s0(:,1); n2=[2:nt1];
dt=ut(end);
tf=dt*N;
% Put u in matrix form, nc rows, N columns
for i=1:N
    u(1:nc,i)=ut((i-1)*nc+1:i*nc);
end
% Forward sequencing and store state histories, s:
for i=1:N
    s(:,i+1)=feval(name,u(:,i),s(:,i),dt,(i-1)*dt,1);
end
% Performance index, terminal constraints & gradients:
[Phi,Phis,Phid]=feval(name,zeros(nc,1),s(:,N1),dt,tf,2);
phi=Phi(1);     % cost function
psi=Phi(n2);    % terminal constraints
la=Phis';       % lecture 14, page 3
% Backward sequencing and store Hu(:,:,i):
for i=N:-1:1
    [fs,fu,fd]=feval(name,u(:,i),s(:,i),dt,(i-1)*dt,3);
    Hu(:,:,i)=la'*fu;
    Phid=Phid+la'*fd;
    la=fs'*la;
end
%---separate into phi and psi derivative terms
phid=Phid(1);
psid=Phid(n2);
%---pull out initial Hu at first time point
Humatrix=Hu(:,:,1);
%---organize into columns (each column is a time step)
for i=2:N,
    Humatrix=[Humatrix Hu(:,:,i)];
end
Humatrix=[Humatrix Phid];
Huphid=Humatrix(1,:)';  %---these are the derivatives of cost function

% (total time) w.r.t. theta (at each time step) and delta_t
% note that derivative w.r.t. delta_t is equal to N

% Hupsid=[Humatrix(n2,:)]'; %---this is the derivative of terminal position
% constraints w.r.t. theta (at each time step) and delta_t
