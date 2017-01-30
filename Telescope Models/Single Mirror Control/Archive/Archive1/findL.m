function L = findL(i,M,p,N_0)
% 

syms x
% f = dot(i,M*i)*x^2 + dot(2*i,(M*p+N_0))*x + dot(p,(M*p+2*N_0));
% f = dot(i,M*i)*x^2 + 2*i.'*(M*p+N_0)*x + dot(p,(M*p+2*N_0));
f = (p+x*i).'*M*(p+x*i) + 2*N_0.'*(p+x*i);
L = solve(f==0,x);
L = double(L);
% L = vpa(L,1000);
L = min(abs(L));
