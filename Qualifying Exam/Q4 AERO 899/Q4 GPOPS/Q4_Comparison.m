%---------------------------------------------------%
% AFIT Qualifying Exam Question #4, 01 Dec 2014     %
% Timothy Coon                                      %
% Advisor: Dr Cobb                                  %
%---------------------------------------------------%
% Q4_Comparison.m
% compare the effects of using different numbers of collocation points and
% different numbers of mesh intervals.

clear all; close all; clc;

N_cols = [2 5 9];
N_mesh = [1 2];

figure('Name','Minimum Energy Path to the Road')
for i = 1:length(N_mesh)
    n_mesh = N_mesh(i);
    for j = 1:length(N_cols)
        n_cols = N_cols(j);
        Q4_Main;
        Q4_PlotComparison
    end
end

