function [ F ] = matrix_concat( A, B, C, varargin )
%matrix_concat combines two or three matrices
%   fx = x-process continuous-time system matrix
%   fy = y-process continuous-time system matrix
%   fz = z-process continuous-time system matrix
%   put the F matrix from each process into one big F matrix. It will be a
%   diagonal matrix with a lot of zeros off-diagonal

if nargin == 2;
    C = [];
end

lx1 = size(A,1);
lx2 = size(A,2);
ly1 = size(B,1);
ly2 = size(B,2);
lz1 = size(C,1);
lz2 = size(C,2);

F = [      A        zeros(lx1,ly2) zeros(lx1,lz2);
     zeros(ly1,lx2)       B        zeros(ly1,lz2);
     zeros(lz1,lx2) zeros(lz1,ly2)       C       ];

end

