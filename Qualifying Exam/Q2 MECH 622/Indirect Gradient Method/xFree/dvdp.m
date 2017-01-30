function [f1,f2,f3]=dvdp(u,s,dt,tf,flg)                   
% Discete Velocity Direction Programming for min tf to specified xf,yf
% s=[x y]'; [th]=control;                 
%
% Modified by Coon 11/2014
%
global c1 c2 xf yf Vmax w Vwx Vwy
x=s(1);
y=s(2);
th = u;
Vriver = calcVriver(y);
dVr_dy = 4*Vmax*(w-2*y)/w^2;
Vwxr = Vwx - Vriver;            % Vel of wind wrt river

if flg==1,     % f1=f(x,th,tf)
    % calculate the next state vector
    f1 = s + [c1*Vwxr*sin(th)+Vriver; c2*Vwy*cos(th)]*dt;
elseif flg==2, % f1=[phi; psi], f2=[phis; psis]', f3=[phid; psid];
    % f1 = [perf index; term constr]; f2 = state grads; f3 = dt grads
    f1=[tf; x-xf; y-yf];
    f2=[0 0; 1 0; 0 1];      % [phi_x phi_y; psix_x psix_y; psiy_x psiy_y] 
    f3=[tf/dt 0 0]';         % [Phi_d]
elseif flg==3, % f1=f_s, f2=f_u, f3=f_dt;
    f1=[1 (-c1*sin(th)+1)*dVr_dy*dt; 0 1];
    f2=dt*[c1*Vwxr*cos(th); -c2*Vwy*sin(th)];
    f3=[c1*Vwxr*sin(th)+Vriver; c2*Vwy*cos(th)];
end

	
