% Q5 P3
% Derive an analytical expression for the diffraction pattern in the focal
% plane caused by a uniform amplitude input (=point source at infinity)
clear; close all; clc;

syms A lm R l lo f1 f2 K u

% fi = [f1 f2];
fi = ones(1,1)*1;
f = 0;                  % initialize surface roughness function

% % roughness term with sines
% for i = 1:length(fi)
%     f = f + fi(i)*sin(pi*i*l/lo);
% end

% roughness term with TSE of sines
O1_tse = 1;
for i = 1:length(fi)
    X = pi*i*l/lo;
    sin_X = 0;
    for k = 0:O1_tse
        sin_X = sin_X + (-1)^k*X^(2*k+1)/factorial(2*k+1);
    end
    f = f + fi(i)*sin_X;
end

th = K*f;

% TSE of exponential
% O2_tse= 1;
% P_tse = 0;
% for k = 0:O2_tse
%     cos_th = (-1)^k*(th^(2*k)/factorial(2*k));
%     sin_th = (-1)^k*(th^(2*k+1)/factorial(2*k+1));
%     P_tse = P_tse + cos_th + 1i*sin_th;
% end

P_tse = exp(1i*K*f);

If = 4*A^2/(lm^2*R^2)*fourier(P_tse,l,u)^2;

f_sub = subs(f,lo,0.1);
If_sub = subs(If,[A,lm,R,lo,K],[1,0.5e-6,0.5,0.1,2*pi/0.5e-6]);

pretty(P_tse)
pretty(If)

figure()
subplot(211)
ezplot(f_sub,[-5e-5,5e-5])
subplot(212)
ezplot(If_sub,[-5e-5,5e-5])
