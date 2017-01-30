function [ d ] = calcD( X, Y, f )
%CALCD uses the surface roughness function to build a discretized surface
%   This is used as the effective path length deviation normal to the
%   gaussian reference sphere (Ref Goodman 6.4.1). In reality, the rough
%   surface will have a more complex effect on the wavefront, but this is a
%   first look.

% assuming a roughness that is radially symmetric, the deviation will be
% defined as a function of the radius from the vertex perpendicular to the
% principle axis. For this, we first need the arc length as a function of
% the radius. (the mirror is a parabola, but the gauss ref sphere
% establishes the normal direction used for the deviation definition here.

% syms lsym

% radius of each point
r = sqrt(X.^2+Y.^2);

% the arclength integral with variable radius was evaluated using mupad.
l0 = sqrt(1+1/f)/2 + f*asinh(1/f)/2;
l = abs(r.*sqrt(1+r.^2/f^2) + f*asinh(r./f))*(1/2);
d = 1e-5*sin(pi*l/l0) + 1e-5*sin(2*pi*l/l0);

% d = 0;              % surface roughness term
% dd = 0;             % derivative of deviation
% di = [1; 1; 1; 1; 1]*1.0e-5;  % weights of deviation (~ wavelengths)
% % build the symbolic function
% for i = 1:length(di)
%     d = d + di(i)*sin(pi*i*lsym/l0);
%     dd = dd + di(i)*cos(pi*i*lsym/l0)*pi*i/l0;
% end
% % evaluate the function
% for j = 1:length(l)
%     D(j,1) = subs(d,lsym,l(j));
%     DD(j,1) = subs(dd,lsym,l(j));
% end
% D = eval(D);
% DD = eval(DD);



end

