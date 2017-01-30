function [Vb,Vbr,Vb_grad] = calcVvector(A,th,Vb0,Vr,varargin)
% CALCVVECTOR calculates the absolute velocity of the boat
% Vb (2x1)= abs velocity of boat (m/s)
% Vbr (2x1) = velocity of boat wrt river
% Vw (2x1)= abs velocity of wind (m/s)
% Vr (2x1)= abs velocity of river (m/s)
% A (1x1)= heading angle of boat wrt +x-axis (rad)
% th (1x1)= angle of the main sail wrt neg boat C/L (rad)
% V0 (2x1)= initial guess for abs velocity (m/s)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% I discovered that, while the accuracy of fsolve() is heavily reliant upon
% the initial guess, fmincon() is far more robust. I made a dummy objective
% file and used an identical function for the constraints to that I used
% with fsolve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global Vw

if nargin == 2  % then we are testing this function
    Aw = deg2rad(130);
    Vw_mag = 10*.514;
    Vw = [Vw_mag*cos(Aw); Vw_mag*sin(Aw)];
    Vr = [0; 0];
%     A = pi/2;
%     th = pi/2;
    Vb0 = [0; 0];
end
c1 = 0.5;     % wind force coeff
c2 = 0.5;     % water drag force coeff

% mu = (c1/c2)^2;                 % force coefficient
rs = [cos((A+pi)+th); sin((A+pi)+th)];    % main sail direction vector

optn = optimset('Display','off');
% [Vb_fsolve,fval] = fsolve(@V_eqns,Vb0,optn);   % does not work as well
% Vb = Vb_fsolve;

[Vb_fmincon,~,~,~,~,Vb_grad_fmincon] = ...
                      fmincon(@dummyJ,Vb0,[],[],[],[],[],[],@V_eqns2,optn);
Vb = Vb_fmincon;
Vb_grad = Vb_grad_fmincon;

Vbr = Vb - Vr;

% Vb = [Vb_fsolve,Vb_fmincon];    % for easy comparison of results


%% fsolve function (not used)
function F = V_eqns(Vb)
    Vwb = Vw - Vb;
    Vbr = Vb - Vr;
    Vbr_hat = Vbr/norm(Vbr);
    sine_alpha = det([rs';Vwb'])/norm(Vwb);
%         sine_theta = det([-Vbr_hat';rs']);
    Fs_mag = c1*norm(Vwb)^2*sine_alpha;
    Fs_br = Fs_mag*sin(th)*Vbr_hat;
    Fd = c2*norm(Vbr)*Vbr;

    F = Fs_br - Fd;
end

%% fmincon NL constraint function
% the equations herein are formed using a force balance with vector
% notation. Satisfying the force balance determines the absolute velocity
% of the boat
function [c,ceq] = V_eqns2(Vb)
    Vwb = Vw - Vb;
    Vbr = Vb - Vr;
%     Vbr_hat = Vbr/norm(Vbr);
    Vbr_hat = [cos(A); sin(A)];
    rs_X_Vwb = det([rs';Vwb']);
    sine_theta = det([-Vbr_hat';rs']);
    Fs_mag = c1*norm(Vwb)*rs_X_Vwb;
    Fs_br = Fs_mag*sine_theta*Vbr_hat;
    Fd = c2*norm(Vbr)*Vbr;

    ceq = Fs_br - Fd;
    c = [];
end

%% fmincon objective function
% this is a dummy function for fmincon() so the constraint equations can be
% solved more robustly than with fsolve()
function J = dummyJ(Vb)
    J = Vb(1);
end
end
