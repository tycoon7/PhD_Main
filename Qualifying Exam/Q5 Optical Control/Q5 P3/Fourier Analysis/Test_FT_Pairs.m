% test fourier transform pairs
% Ref Gaskill Table 7-7
clear; close all; clc;

syms x xi

% exp((+/-)*1i*pi*(x^2-(1/8))) => exp((-/+)*1i*pi*(xi^2-(1/8)))

f1 = exp(1i*pi*(x^2-(1/8)));
F1_gaskill = exp(-1i*pi*(xi^2-(1/8)));
F1 = fourier(f1,x,xi);

pretty(F1_gaskill)
pretty(F1)
