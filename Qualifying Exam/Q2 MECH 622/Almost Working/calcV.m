function V = calcV(u,Vriver,chi0)
% calcV holds the set of equations to calculate relative velocity of the
% boat and other stuff using fsolve
% u = inputs
% Vriver = velocity of the river (x-direction)
% chi0 = initial guess for chi = [v Wr alpha]

global Vwx Vwy

th = u(1);
A = u(2);
Aw = atan(Vwy/(Vwx-Vriver));
psi = Aw - A;
W = sqrt((Vwx-Vriver)^2 + Vwy^2);

% guess a value of mu (will depend on boat properties)
mu = 1;         % coefficients

options = optimset('Display','off');
[X fval] = fsolve(@V_Wr_alpha_eqns,chi0,options);
V = X(1);
 
%% nested function

    function F = V_Wr_alpha_eqns(x)
        v = x(1);
        Wr = x(2);
        alpha = x(3);

        F = [mu^2*Wr^2*sin(alpha)*sin(th) - v^2;
             W*sin(psi) - Wr*sin(alpha+th);
             v^2 + W^2 - 2*v*W*cos(psi) - Wr^2];
    end
end
 

