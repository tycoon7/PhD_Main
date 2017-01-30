%---------------------------------------------------%
% AFIT Qualifying Exam Question #2, 02 Dec 2014     %
% Timothy Coon                                      %
% Advisor: Dr Cobb                                  %
%---------------------------------------------------%
% Q2_Main.m                                         %
% The problem solved here is given as follows:      %
% Minimize the time to get across the river in a    %
% sailboat                                          %
%    x_dot = vx                                     %
%    y_dot = vy                                     %
% and the boundary conditions                       %
%    x0 = 85, xf = 115                              %
%    y0 = 0, yf = 200                               %
%---------------------------------------------------%
global Vmax w Vw

Vmax = 32*0.514;            % (m/s) max river velocity
w = 200;                    % (m) river width
Aw = deg2rad(130);          % (rad) angle of the wind wrt +x-axis
Vw_mag = 10*.514;           % (m/s) magnitude of wind velocity
Vw = [Vw_mag*cos(Aw); Vw_mag*sin(Aw)];  % (m/s) wind veloctity vector

% tune boundary conditions
t0 = 0;                         % (m) force start at time = zero
tfmin = 100; tfmax = 300;

% state limits
xmin = 0; xmax = 500;           % (m) min/max y-position
ymin = 0; ymax = 200;           % (m) min/max y-position
% control limits
Tmin = -pi/2; Tmax = pi/2;      % (rad) min/max sail angle wrt -C/L
Amin = -2*pi; Amax = 2*pi;      % (rad) min/max heading angle
% state boundary conditions
x0 = 85; xf = 115;                % (m)
y0 = 0;  yf = 200;


%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
% time
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tfmin; 
bounds.phase.finaltime.upper = tfmax;

% states
% fixed initial state
bounds.phase.initialstate.lower = [x0 y0]; 
bounds.phase.initialstate.upper = [x0 y0];
% bounds for states after initial
bounds.phase.state.lower = [xmin ymin]; 
bounds.phase.state.upper = [xmax ymax];
% free final state
bounds.phase.finalstate.lower = [xf yf]; 
bounds.phase.finalstate.upper = [xf yf];

% control limits
bounds.phase.control.lower = [Tmin Amin]; 
bounds.phase.control.upper = [Tmax Amax];

% path constraints
% bounds.phase.path.lower = [0 0];
% bounds.phase.path.upper = [0 0];
% The integral constraints are inequalities. Express this by setting the
% lower bound to be the minimum required, then a guess for the upper
% bound which is NOT ARBITRARILY HIGH, MUST BE "CLOSE". 
% (s) is the number of UGSs and, hence, the number of integrals.
% bounds.phase.integral.lower = [100];
% bounds.phase.integral.upper = [300];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t0; (tfmax+tfmin)/2];
guess.phase.state   = [[x0; xf],[y0; yf]];
guess.phase.control = [[-pi/4; -pi/4],[pi/6; pi/6]];
% guess.phase.integral = [120];

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name = 'Q2-Problem';
setup.functions.continuous = @Q2_Continuous;
setup.functions.endpoint = @Q2_Endpoint;
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparseNaN';
setup.scales.method = 'none';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-2;
setup.mesh.maxiterations = 10;
n_mesh = 10; n_cols = 10;
setup.mesh.phase.fraction = (1/n_mesh)*ones(1,n_mesh);
setup.mesh.phase.colpoints = n_cols*ones(1,n_mesh);
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
setup.nlp.ipoptoptions.tolerance = 1e-2;
setup.nlp.ipoptoptions.maxiterations = 50;
setup.displaylevel = 2;

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
output.result.nlptime
solution = output.result.solution;
GPOPS2_ExitFlags(output)
