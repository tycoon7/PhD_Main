classdef MirrorClass < handle
    %MIRRORCLASS for mirror objects
    %   Detailed explanation goes here

%%
    properties
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
        
    end
        
            
end

