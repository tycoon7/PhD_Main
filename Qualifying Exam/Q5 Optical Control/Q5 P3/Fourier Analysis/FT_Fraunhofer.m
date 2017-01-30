% fourier tranform

syms x u k R c1

f = exp(1i*k*c1*sin(pi*x/R));
fourier(f,x,u)