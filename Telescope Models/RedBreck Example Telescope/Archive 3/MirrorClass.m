classdef MirrorClass < handle
    %MIRRORCLASS for mirror objects
    %   Detailed explanation goes here

%%
    properties
        % raytrace properties
%         P;      % start point, base coords
%         I;      % incident ray direction, base coords
        p;      % incident ray start point, mirror coords
        i;      % incident ray direction, mirror coords
        r;      % reflected ray direction, mirror coords
        rho;    % intersection point, mirror coords

        % mirror parameters
        e;      % mirror surface eccentricity
        f;      % mirror focal length
        V;      % mirror vertex position, base coords
        R_mb;   % mirror frame rotation matrix
        M;      % surface dyadic
        psi;    % principle axis direction vector
        N0;     % radius of curvature vector
    end

%%
    methods
        function obj = MirrorClass(e,f,V,psi,R_mb)
            if nargin > 0
                obj.e = e;
                obj.f = f;
                obj.V = V;
                obj.R_mb = R_mb;
                obj.psi = psi;
                obj.M = obj.findM(obj.e,obj.psi);
                obj.N0 = obj.findN0(obj.f,obj.e,obj.psi);
            end
        end        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [P2,I2] = reflect(obj,P1,I1)
            obj.p = obj.R_mb*(P1-obj.V);
            obj.i = obj.R_mb*I1;
            obj.L = obj.findL(obj.i,obj.p,obj.N0,obj.M,obj.psi);
            obj.rho = obj.findRho(obj.p,obj.L,obj.i);
            obj.r = obj.findr(obj.f,obj.psi,obj.i,obj.M,obj.N0,obj.rho);
            P2 = obj.R_mb'*obj.rho + obj.V;
            I2 = obj.R_mb'*obj.r;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function r = findr(obj,f,psi,i,M,N0,rho)
            if f ~= inf
                obj.N_vec = N0 + M*rho;       % surface normal vector at intersection
                obj.N_hat = -sign(dot(i,obj.N_vec))*obj.N_vec/norm(obj.N_vec);
            else
                obj.N_hat = psi;
            end
            R = eye(3)-2*obj.N_hat*transpose(obj.N_hat);
            r = R*i; 
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static, Access = private)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function M = findM(e,psi)
            M = eye(3)-e^2*(psi*psi');
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function N0 = findN0(f,e,psi)
            if f == inf
                N0 = psi;
            else
                N0 = -f*(1+e)*psi;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function L = findL(i,p,N0,M,psi)
            syms x
            if ~isequal(N0,psi) 
                F = (p+x*i).'*M*(p+x*i) + 2*N0.'*(p+x*i);
            else
                F = psi'*(p+i*x);
            end
            L = solve(F==0,x);
            L = double(L);
            % L = vpa(L,1000);
            L = min(abs(L));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = findRho(p,L,i)
            rho = p+L*i;
        end
        
    end
        
            
end

