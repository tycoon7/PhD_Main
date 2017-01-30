%---------------------------------------------------%
% Single Mirror Control:                            %
% AFIT AERO 899 Summer 2014 Dr Cobb                 %
% Timothy Coon                                      %
%---------------------------------------------------%
% The problem solved here is given as follows:      %
% (see Burl Ch 2, C.E. 2.1)                         %
%   Minimize spot size                              %
% subject to the dynamic constraints                %
%    dx1/dt = x2                                    %
%    dx2/dt = (1/J)*(-G*x2 + u + w)                 %   
%    x1 = angle error                               %
%    x2 = angle error rate                          %
% and the boundary conditions                       %
%    x1(0) = 0, x2(0) = 0                           %
%    x1(t_f) = FREE, x2(t_f) = FREE                 %
%---------------------------------------------------%

clear all;
close all; clc;

auxdata.f = 0.70;              % (m) focal length of parabolic mirror
auxdata.G = 0.05;             % (kg-m) friction force per (m/s)
auxdata.J = 4*0.2/(3*pi);   % (kg-m^2) MOI of primary mirror
% auxdata.w = 3.*randn(100,1);     % (N) random disturbance force
% auxdata.ps = pi*rand;         % (rad) random phase shift
auxdata.ps = 0;
                                
% time boundary conditions
t0 = 0;                       % (s) force start at time = zero
tf = 1;                       % (s) fixed simulation time
% auxdata.tw = linspace(t0,tf,length(auxdata.w));

% state limits
x1min = -2; x1max = 2;        % (m) min/max angle error
x2min = -5; x2max = 5;        % (m) min/max angle error rate
fumin = -5; fumax = 5;        % (N) min/max control force

% control limits
% umin = -1; umax = 1;          % (N-m) control torque
umin = -9; umax = 9;          % (W) control power

% state boundary conditions
x10 = 0;                      % (m)
x20 = 0;                      % (m)
fu0 = 0;                      % (N)

%-------------------------------------------------------------------------%
%----------------------- Setup for Problem Bounds ------------------------%
%-------------------------------------------------------------------------%
iphase = 1;
% time
bounds.phase.initialtime.lower = t0; 
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf; 
bounds.phase.finaltime.upper = tf;

% states
% fixed initial state
bounds.phase.initialstate.lower = [x10,x20,fu0]; 
bounds.phase.initialstate.upper = [x10,x20,fu0];
% bounds for states after initial
bounds.phase.state.lower = [x1min,x2min,fumin]; 
bounds.phase.state.upper = [x2max,x2max,fumax];
% free final state
bounds.phase.finalstate.lower = [x1min,x2min,fumin]; 
bounds.phase.finalstate.upper = [x1max,x2max,fumax];

% control limits
bounds.phase.control.lower = [umin]; 
bounds.phase.control.upper = [umax];

% path constraints
% bounds.phase.path.lower = 1;
% bounds.phase.path.upper = 1;
% The integral constraints are inequalities. Express this by setting the
% lower bound to be the minimum required, then a guess for the upper
% bound which is NOT ARBITRARILY HIGH, MUST BE "CLOSE". 
bounds.phase.integral.lower = [0];
bounds.phase.integral.upper = [100];

%-------------------------------------------------------------------------%
%---------------------- Provide Guess of Solution ------------------------%
%-------------------------------------------------------------------------%
guess.phase.time    = [t0; tf];
guess.phase.state   = [[x10, x10, fu0];[x20, x20, fu0]];
guess.phase.control = [[0; 0]];
% the integral guess is to be a row vector (1xnd)
guess.phase.integral = [0];

%-------------------------------------------------------------------------%
%------------- Assemble Information into Problem Structure ---------------%        
%-------------------------------------------------------------------------%
setup.name = 'SMC-Problem';
setup.functions.continuous = @SMC_Continuous;
setup.functions.endpoint = @SMC_Endpoint;
setup.auxdata = auxdata;
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

setup.mesh.phase.fraction = 0.05*ones(1,20);
setup.mesh.phase.colpoints = 5*ones(1,20);

setup.nlp.solver = 'ipopt';
setup.nlp.ipoptoptions.linear_solver = 'mumps';
setup.nlp.ipoptoptions.tolerance = 1e-4;
setup.nlp.ipoptoptions.maxiterations = 3000;
% setup.nlp.snoptoptions.tolerance = 1e-6;
% setup.nlp.snoptoptions.maxiterations = 2000;

setup.displaylevel = 2;

setup.mesh.colpointsmin = 10;
setup.mesh.colpointsmax = 10;
% setup.method = 'RPM-Integration';
setup.mesh.phase.colpoints = 5*ones(1,20);

% setup.scales.method = 'automatic-bounds';

%-------------------------------------------------------------------------%
%------------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------------%
output = gpops2(setup);
output.result.nlptime
solution = output.result.solution;
GPOPS2_ExitFlags(output)

%% -----------------------------------------------------------------------%
%------------------------------- Plot Solution ---------------------------%
%-------------------------------------------------------------------------%
SMC_Plot
