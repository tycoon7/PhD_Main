function [f1,f2,f3]=dvdp(u,s,dt,t,flg)                   
% Discete Velocity Direction Programming for min tf to specified xf,yf

% s=[u v y x]'; [th A]=control;                 
%
% Modified by Coon 11/2014
%
global xf yf Vmax Vwx Vwy
x=s(1);
y=s(2);
th = u(1);
A = u(2);
co=cos(A);
si=sin(A);
Vriver = calcVriver(y);
% options = optimset('Jacobian','off');
options = optimoptions('fsolve','Jacobian','off');
X0 = [0; 0; 0];     % [V; Wr; alpha]
[X,fval,exitflag,output,J] = fsolve(@(X)calcV(X,u,Vriver),X0,options);
% J = [ dV/dV  dV/dWr  dV/da;
%      dWr/dV dWr/dWr dWr/da;
%       da/dV  da/dwr  da/da];
V = X(1);       % relative velocity of boat wrt river

if flg==1,     % f1=f(x,th,tf)
    % calculate the next state vector
    f1=s+dt*[V*co+Vriver; V*si;];
elseif flg==2, % f1=[phi; psi], f2=[phis; psis], f3=[phid; psid];
    % f1 = [perf index; term constr]; f2 = state grads; f3 = dt grads
    f1=[t; x-xf; y-yf];
    f2=[0 1 0; 0 0 1];      % [Phi_x Phi_y]
    f3=[t/dt 0 0]';         % [Phi_d]
elseif flg==3, % f1=f_s, f2=f_u, f3=f_dt;
    f1=[1 0; 0 0];
    f2=dt*[-si; -dt*si/2; co; dt*co/2];
    f3=[co; u+dt*co; si; V+dt*si];
end

	
