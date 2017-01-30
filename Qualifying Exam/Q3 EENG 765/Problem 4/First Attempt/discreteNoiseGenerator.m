function [ wk ] = discreteNoiseGenerator( Qk, num_sources, num_steps, G )
%DISCRETENOISEGENERATOR Generates discrete noise from Qd using Cholesky
%Factorization method
%   Qk = discrete covariance matrix
%   num_sources = number of white noise sources
%   num_steps = number of time steps in the process
% 
%   Reference Brown and Hwang, Section 3.10 for details. The goal is to
%   generate (m) vectors of white gaussian sequences to be used as discrete
%   noise in simulation, where (m) is the number of input noise sources. To
%   generate these vectors, begin with vectors of normally distributed
%   unitary noise sequence vectors (uk) and operate on them with a linear
%   operator, [Ck], chosen to yield (wk) with the desired covariance
%   structure found in [Qk]. [Ck] is not unique, but a simple way to
%   generate a suitable [Ck} is to assume lower triangular and use Cholesky
%   Factorization. Built-in MATLAB function chol() makes this painless even
%   for large matrices.
%
%   Added: test if the input Qk is singular and use singular value
%   decomposition method if it is.

[m n] = size(Qk);
R = rank(Qk);
test_singular = R < min(m,n);
if test_singular == false
    Ck = chol(Qk,'lower');      % Ck is lower triangular linear operator
else
    [U T V] = svd(Qk);
    S = sqrt(T);
    Ck = U*S;                   % Ck is not lower triangular
end
    
% this is what I had before, but I didn't take into account that G is used
% to generate Qk....
% uk = G*randn(num_sources,num_steps);  % normally distributed unitary sequences

uk = randn(num_sources,num_steps);  % normally distributed unitary sequences
G_sign = G./G;
G_sign(isnan(G_sign)) = 0;
wk = Ck*G_sign*uk;

% check that the C matrix satisfies Qd = CC' = cov(wk') (gets close)
Qk_check = cov(wk');    % definition of Qd
check = Qk-Qk_check;    % could use this for error-checking

end

