classdef RayClass < handle
    %RAYCLASS for ray objects
    %   Detailed explanation goes here
    
    properties
        P;          % ray start point, base coords
        I;          % ray direction, base coords
        rayNum;     % ray number
        seg = struct('P',{},'I',{},'p',{},'i',{},'r',{},'N_vec',{},'N_hat',{},'L',{});
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = RayClass(obj,rayNum,P0,I0,mirObj)
            % constructor 
            % P0 = start point of chief ray
            % I0 = start direction of all rays
            obj.rayNum = rayNum;
            obj.P = RayClass.findP(rayNum,P0);
            obj.I = I0;
            obj.rayTrace(mirObj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rayTrace(obj,mirObj)
            for j = 1:length(obj,mirObj)
                obj.seg(j).p = mirObj(j).R_mb*(obj.seg(j).P-mirObj(j).V);
                obj.seg(j).i = mirObj(j).R_mb*obj.seg(j).I;
                obj.seg(j).L = mirObj(j).findL(obj.seg(j).p,obj.seg(j).i,mirObj(j).N0,mirObj(j).M,mirobj(j).psi);
                obj.seg(j).rho = obj.findRho(obj.seg(j).p,obj.seg(j).i,obj.seg(j).L);
                obj.seg(j) = obj.findNRr(obj.seg(j),mirObj(j));
                obj.seg(j+1).P = mirObj(j).R_mb'*obj.seg(j).rho + mirObj(j).V;
                obj.seg(j+1).I = mirObj(j).R_mb'*obj.seg(j).r;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function seg = findNRr(obj,seg,mirObj)
            if mirObj.f ~= inf
                seg.N_vec = mirObj.N0 + mirObj.M + seg.rho;
                seg.N_hat = -sign(dot(obj.i,obj.N_vec))*obj.N_vec/norm(obj.N_vec);
            else
                seg.N_hat = psi;
            end
            R = eye(3)-2*seg.N_hat*transpose(seg.N_hat);
            seg.r = R*seg.i; 
        end
    end
    
    methods (Access = private, Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function P = findP(rayNum,P0)
            % find startpoint of ray rayNum
            % make this work for any number of rays and specify beamwidth
            switch rayNum
                case 1
                    P = P0;
                case 2
                    P = P0 + [0.00;  0.05;  0.00];
                case 3
                    P = P0 + [0.00; -0.05;  0.00];
                case 4
                    P = P0 + [0.00;  0.00;  0.05];
                case 5
                    P = P0 + [0.00;  0.00; -0.05];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function L = findL(p,i,N0,M,psi)
            syms x
            if ~isequal(N0,psi) 
                F = (p+x*i).'*M*(p+x*i) + 2*N0.'*(p+x*i);
            else
                F = psi'*(p+i*x);
            end
            L = solve(F==0,x);
            L = double(L);
            L = min(abs(L));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rho = findRho(p,i,L)
            rho = p+L*i;
        end
        
    end
end

