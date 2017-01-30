classdef BeamClass < handle
    %BEAMCLASS draws rays to match beam specifications
    %   The beam has a specified width and ray density. BeamClass uses this
    %   to find the grid of ray start points (P0's) and assign it to a ray
    %   in an indexing fashion. The line-of-sight matrices, C_LOS, for each
    %   ray are combined here to produce the output wavefront matrix, C_WF
    
    %   I could make a lot of beam objects at different field angles with
    %   different weights and build any input I want to (angular spectrum)
    
    properties
        BW;                     % user-specified beamwidth
        rayDensity;             % user-specified density of the rays
        P0_chief;               % start point of chief ray
        I0;                     % direction vector of chief ray
        mirObj;                 % mirror object
        nrays;                  % number of rays
        nmirs;                  % number of mirrors
        nsegs;                  % number of ray segments
        ray@RayClass;           % ray objects (Kool!!)
        C_WF;                   % wavefront matrix
        
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = BeamClass(BW,rayDensity,P0_chief,I0,mirObj)
            obj.BW = BW;
            obj.rayDensity = rayDensity;
            obj.P0_chief = P0_chief;
            obj.I0 = I0;
            obj.mirObj = mirObj;
            obj.nmirs = length(mirObj);
            obj.nsegs = obj.nmirs;
            rayStartGrid = obj.calcRayStartPts(BW,rayDensity,I0);
            P0 = bsxfun(@plus,P0_chief,rayStartGrid);
            obj.nrays = length(rayStartGrid);
            for r = 1:obj.nrays
                obj.ray(r) = RayClass(r,P0(:,r),I0,mirObj);
            end
%             obj.calcC_WF;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = calcC_WF(obj)
            obj.C_WF = obj.ray(1).C_los;
            for i = 2:obj.nrays
                obj.C_WF = [obj.C_WF; obj.ray(i).C_los];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods (Access = public, Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotNomRaytrace(obj,fhandle)
            figure(fhandle)
            hold on
            for s = 1:obj.nsegs
                for r = 1:obj.nrays
                    P(r,1) = obj.ray(r).seg(s).P(1);
                    P(r,2) = obj.ray(r).seg(s).P(2);
                    P(r,3) = obj.ray(r).seg(s).P(3);
                    % intersection points
%                     plot3(obj.ray(r).seg(s).P(1),obj.ray(r).seg(s).P(2),obj.ray(r).seg(s).P(3),'g+','MarkerSize', 10);
                    plot3(P(r,1),P(r,2),P(r,3),'g+','MarkerSize', 10);
                    % ray direction
                    if s ~= obj.nsegs
%                     quiver3(obj.ray(r).seg(s).P(1),obj.ray(r).seg(s).P(2),obj.ray(r).seg(s).P(3),...
%                             obj.ray(r).seg(s).I(1),obj.ray(r).seg(s).I(2),obj.ray(r).seg(s).I(3),'Color','r')
                    arrow([obj.ray(r).seg(s).P(1),obj.ray(r).seg(s).P(2),obj.ray(r).seg(s).P(3)],...
                          [obj.ray(r).seg(s+1).P(1),obj.ray(r).seg(s+1).P(2),obj.ray(r).seg(s+1).P(3)],...
                           'EdgeColor','r','FaceColor','r')
                    arrow([P(r,1),P(r,2),P(r,3)],...
                          [obj.ray(r).seg(s+1).P(1),obj.ray(r).seg(s+1).P(2),obj.ray(r).seg(s+1).P(3)],...
                           'EdgeColor','r','FaceColor','r')
                    end
                end
%                 griddedInterpolant(P(:,1),P(:,2),P(:,3))
            end
            for j=1:obj.nmirs
                plot3(obj.mirObj(j).V(1),obj.mirObj(j).V(2),obj.mirObj(j).V(3),'.','MarkerSize',30);
            end
            hold off
            axis equal
            xlabel('x'); ylabel('y');
            view(3)
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotNomSpots(obj,fhandle)
            figure(fhandle)
            title('Spot Diagram')
            hold on
            % Plot spot diagram in the structural xz-plane
            parfor nr = 1:obj.nrays
                plot(obj.ray(nr).seg(end).P(1),obj.ray(nr).seg(end).P(3),'b.','MarkerSize',10);
            end
            hold off
%             grid(ax,'on')
            axis equal
            xlabel('x'); ylabel('z');
%             legend('Exact')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotPertSpots(obj,fhandle,gamma)
            figure(fhandle)
            title('Spot Diagram')
            errScale = 1;
            gamma = errScale*gamma;
            hold on
            % Plot spot diagram in the structural xz-plane
            parfor nr = 1:obj.nrays
                plot(obj.ray(nr).seg(end).P(1)-gamma(1,nr),obj.ray(nr).seg(end).P(3)-gamma(3,nr),'r+','MarkerSize',10);
            end
            grid on
            hold off
            xlabel('x'); ylabel('z');
%             legend('Nominal','Perturbed')
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function plotDiffractionPattern(obj,fhandle,integralType)
            figure(fhandle)
            % build OPD "function" on Gaussian Reference Sphere
            % I want to do this both as a discrete FFT and analytically
            % with a polynomial curve fit (probably in two diff methods)
            switch integralType   % maybe I only need Fraunhofer??
                case 'Rayleigh-Sommerfeld'
                    % stuff
                case 'Fresnel'
                    % stuff
                case 'Fraunhofer'
                    % stuff
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
    methods (Access = private, Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function rayStartGrid = calcRayStartPts(BW,rayDensity,I0)
            % find the coordinates of the points in a circular grid
            % make this in XY-plane and use a rotation angle to normal with
            % input ray direction, I0
            R = BW/2;
            rd1 = rayDensity/2;     % radial ray density
            rd2 = 3*rayDensity;     % circular ray density
            if numel(R) > 1;
                error('R?')
            end
            if numel(rd1) > 1;
                error('rd1?')
            end
            [r,phi] = meshgrid(linspace(R/rd1,R,rd1),...
                               linspace(2*pi/rd2,2*pi,rd2));
            if I0 == [1; 0; 0]
            % make this more dynamic to work with any I0 direction
                X = zeros(size(r));
                Y = r.*cos(phi);
                Z = r.*sin(phi);
            elseif I0 == [0; -1; 0]
                X = r.*sin(phi);
                Y = zeros(size(r));
                Z = r.*cos(phi);
            else
                error('Input ray direction not recognized, unable to make field')
            end
            origin = [0; 0; 0];
            rayStartGrid = [origin, [X(:)';Y(:)';Z(:)']];
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    
end

