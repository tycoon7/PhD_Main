function [ t,theta,u,w ] = mirrorDynamics( )
%MIRRORDYNAMICS Summary of this function goes here
% Burl Ch 2 Computer Exercise 2.1
% Altitude angle error antenna
% clear all; close all; clc;

% time vector
dt = 0.02;
tf = 0.5;
t = 0:dt:tf;

% parameters
G = 5;      % kg-m^2/s
J = 100;   % kg-m^2
x0 = [0 0];

% system matrices
A = [0 1; -20/J -(80+G)/J];      % closed-loop system matrix
B = [0; 1];
C = [-20 -80];
D = [];

% disturbance input
w = 10*randn(length(t),1);
w = w-mean(w);              % remove bias

% linear system
sys = ss(A,B,C,D);

% Simulation
[y,t,x] = lsim(sys,w,t,x0);
% u = -2000*x(:,1)-800*x(:,2);
% lsim(sys,w,t,x0)
theta = x(:,1);
u = y;
end

