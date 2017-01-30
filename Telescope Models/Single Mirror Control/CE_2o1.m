% Burl Ch 2 Computer Exercise 2.1
% Altitude angle error antenna
clear all; close all; clc;

% time vector
t = 0:0.01:10;

% parameters
G = 5;      % kg-m^2/s
J = 1000;   % kg-m^2
x0 = [0 0];

% system matrices
A = [0 1; -2000/J -(800+G)/J];      % closed-loop system matrix
B = [0; 1];
C = [1 0; 0 1];
D = [];

% disturbance input
w = randn(length(t),1);

% linear system
sys = ss(A,B,C,D);

% Simulation
lsim(sys,w,t,x0)
[y,t,x] = lsim(sys,w,t,x0)