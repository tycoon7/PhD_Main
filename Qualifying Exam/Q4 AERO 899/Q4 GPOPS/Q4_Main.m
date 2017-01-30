%---------------------------------------------------%
% AFIT Qualifying Exam Question #4, 01 Dec 2014     %
% Timothy Coon                                      %
% Advisor: Dr Cobb                                  %
%---------------------------------------------------%
% Q4_Main.m                                         %
% The problem solved here is given as follows:      %
% Minimize effort to get to the road subject to the %
% dynamic constraints (x is ind. var.)              %
%    dy = (dy/dx)*dx                                % 
% and the boundary conditions                       %
%    x0 = 0, y(0) = 1                               %
%    xf = 3, y(xf) = FREE                           %
%---------------------------------------------------%

% x boundary conditions
x0 = 0;                       % (l) force start at time = zero
xf = 3;                       % (l) fixed simulation time
% state limits
ymin = 0; ymax = 3;         % (l) min/max y-position
% control limits
umin = -10; umax = 10;        % (-) min/max slope
% state boundary conditions
y0 = 1;                      % (l)

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
% time = x-position
bounds.phase.initialtime.lower = x0; 
bounds.phase.initialtime.upper = x0;
bounds.phase.finaltime.lower = xf; 
bounds.phase.finaltime.upper = xf;

% states
% fixed initial state
bounds.phase.initialstate.lower = [y0]; 
bounds.phase.initialstate.upper = [y0];
% bounds for states after initial
bounds.phase.state.lower = [ymin]; 
bounds.phase.state.upper = [ymax];
% free final state
bounds.phase.finalstate.lower = [ymin]; 
bounds.phase.finalstate.upper = [ymax];

% control limits
bounds.phase.control.lower = [umin]; 
bounds.phase.control.upper = [umax];

% path constraints
bounds.phase.integral.lower = [0];
bounds.phase.integral.upper = [100];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [x0; xf];
guess.phase.state   = [y0;y0];
guess.phase.control = [0; 0];
guess.phase.integral = [0];

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name = 'Q4-Problem';
setup.functions.continuous = @Q4_Continuous;
setup.functions.endpoint = @Q4_Endpoint;
setup.bounds = bounds;
setup.guess = guess;
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.derivatives.dependencies = 'sparseNaN';
setup.scales.method = 'none';
setup.method = 'RPM-Differentiation';
setup.mesh.method = 'hp-PattersonRao';
setup.mesh.tolerance = 1e-4;
setup.mesh.maxiteration = 10;
setup.mesh.phase.fraction = (1/n_mesh)*ones(1,n_mesh);
setup.mesh.phase.colpoints = n_cols*ones(1,n_mesh);
setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
setup.nlp.ipoptoptions.tolerance = 1e-4;
setup.nlp.ipoptoptions.maxiterations = 100;
setup.displaylevel = 2;

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
output.result.nlptime
solution = output.result.solution;
GPOPS2_ExitFlags(output)
