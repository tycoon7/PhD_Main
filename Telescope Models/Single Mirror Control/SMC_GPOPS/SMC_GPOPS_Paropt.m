% SMC_Paropt
% This script runs SMC_GPOPS code multiple times with varying focal lengths
% and saves the relevant data

% clear all;
% close all; clc;

fLength = 0.65:0.01:0.75;       % focal length parameters
SMC_Setup;                      % setup GPOPS problem

for i = 1:length(fLength)
    auxdata.f = fLength(i);
    output(i) = gpops2(setup);
    ExitFlag(i) = GPOPS2_ExitFlags(output(i));
    solution = output(i).result.solution;
    SMC_Plot
end