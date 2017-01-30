classdef RayClass < handle
    %RAYCLASS for ray objects
    %   Detailed explanation goes here
    
    properties
        rayNum;     % ray number
        Clos;       % ray output matrix
        seg = struct('P',{},'I',{},'p',{},'i',{},'r',{},'N_vec',{},...
                     'N_hat',{},'L',{},'R',{},'dwdw',{},'dwdu',{});
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = RayClass(rayNum,P0,I0,mirObj)
            % constructor 
            % P0 = start point of chief ray
            % I0 = start direction of all rays
            obj.rayNum = rayNum;
            obj.seg(1).P = RayClass.findP(rayNum,P0,I0);
            obj.seg(1).I = I0;
            obj.rayTrace(mirObj);
%             obj.Clos = outMatrix(obj);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rayTrace(obj,mirObj)
            for j = 1:length(mirObj)
                obj.seg(j).p = mirObj(j).R_mb*(obj.seg(j).P-mirObj(j).V);
                obj.seg(j).i = mirObj(j).R_mb*obj.seg(j).I;
                obj.seg(j).L = obj.findL(obj.seg(j).p,obj.seg(j).i,mirObj(j).N0,mirObj(j).M,mirObj(j).psi);
                obj.seg(j).rho = obj.findRho(obj.seg(j).p,obj.seg(j).i,obj.seg(j).L);
                obj.seg(j) = obj.findNRr(obj.seg(j),mirObj(j));
%                 obj.seg(j) = obj.sensitivities(obj.seg(j),mirObj(j));
                obj.seg(j+1).P = mirObj(j).R_mb'*obj.seg(j).rho + mirObj(j).V;
                obj.seg(j+1).I = mirObj(j).R_mb'*obj.seg(j).r;
            end
        end
        
    end
    
    methods (Access = private, Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function P = findP(rayNum,P0,I0)
            % find startpoint of ray rayNum
            % make this work for any number of rays and specify beamwidth
            switch rayNum
                case 1
                    P = P0;
                case 2
                    P = P0 + [ 0.05;  0.00;  0.00];
                case 3
                    P = P0 + [-0.05;  0.00;  0.00];
                case 4
                    P = P0 + [ 0.00;  0.00;  0.05];
                case 5
                    P = P0 + [ 0.00;  0.00; -0.05];
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function seg = findNRr(seg,mirObj)
            if mirObj.f ~= inf
                seg.N_vec = mirObj.N0 + mirObj.M*seg.rho;
                seg.N_hat = -sign(dot(seg.i,seg.N_vec))*seg.N_vec/norm(seg.N_vec);
            else
                seg.N_hat = mirObj.psi;
                seg.N_vec = mirObj.psi;
            end
            seg.R = eye(3)-2*seg.N_hat*transpose(seg.N_hat);
            seg.r = seg.R*seg.i; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function seg = sensitivities(seg,mirObj)
            dr_di = drdi(seg.N_hat,seg.i,seg.N_vec,mirObj.M,seg.L,seg.R);
            dr_da = drda(seg.N_hat,seg.i,seg.N_vec,mirObj.M);
            dr_dL = zeros(3,1);
            dr_dt = drdt(mirObj.e,seg.i,seg.r,seg.N_hat,seg.N_vec,mirObj.M,seg.p,mirObj.q);
            dr_dd = drdd(seg.i,seg.N_hat,seg.N_vec,mirObj.M);
            
            dg_di = dgdi(seg.L,seg.R);
            dg_da = dgda(seg.R);
            dg_dL = zeros(3,1);
            dg_dd = dgdd(seg.r,seg.i,seg.N_hat);
            dg_dt = dgdt(seg.r,seg.i,seg.N_hat,seg.p,mirObj.q);
            
            dL_di = zeros(1,3);
            dL_da = zeros(1,3);
            dL_dL = 1;
            dL_dd = dLdd(seg.r,seg.i,seg.N_hat);
            dL_dt = dLdt(seg.i,seg.r,seg.N_hat,seg.p,mirObj.q);
            
            seg.dwdw = [dr_di dr_da dr_dL; dg_di dg_da dg_dL; dL_di dL_da dL_dL];
            seg.dwdu = [dr_dt dr_dd; dg_dt dg_dd; dL_dt dL_dd];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         function ray = outMatrix(ray)
%             ray.Clos = 
    end
end

