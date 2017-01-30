% Qualifying Exam - EENG 765
% Problem 3b
% Use the Van Loan Method to calculate the covariace matrix and determine
% if state variable one and state variable two are uncorrelated

syms s1 s2 b

%% Given
% system matrix
F = [1 0 1; 
     0 1 1; 
     1 0 0];
% PSD of noise
W = eye(2);
% coefficient matrix for unity white noise disturbances
G = [2*sqrt(2*s1^2*b)        0       ;
            0          sqrt(2*s2^2*b); 
            0                0       ];
% assume step size of one
dt = 1;

%% Van Loan Method
A = [F,        G*W*G.'; 
     zeros(3),   F.'  ];
B = expm(A);
phi = B(1:3,1:3);
Q = phi*B(1:3,4:6);
s1 = 1; s2 = 2; b = 1;
Q_num = eval(subs(Q));
% replace all nonzero entries with ones
Qz = spones(Q);

%% Conclusion
% There are no nonzero entries in the covariance matrix, thus all states
% are somewhat correlated because the processes are assumed zero-mean