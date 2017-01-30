function [ t,theta ] = OL_mirrorDynamics( )
%MIRRORDYNAMICS Summary of this function goes here
% Burl Ch 2 Computer Exercise 2.1
% Altitude angle error antenna
% clear all; close all; clc;

% time vector
t = 0:0.02:0.5;

% parameters
G = 5;      % kg-m^2/s
J = 1000;   % kg-m^2
x0 = [0 0];

% system matrices
A = [0 1; 0 -G/J];      % closed-loop system matrix
B = [0; 1];
C = [1 1];
D = [];

% disturbance input
w = randn(length(t),1);
w = w-mean(w);

% linear system
sys = ss(A,B,C,D);

% Simulation
[y,t,x] = lsim(sys,w,t,x0);
lsim(sys,w,t,x0)
theta = y(:,1);
end

