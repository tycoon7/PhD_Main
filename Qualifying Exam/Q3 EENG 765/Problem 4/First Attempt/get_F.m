function [ F ] = matrix_concat( fx, fy, fz, varargin )
%GET_F builds the big F matrix
%   fx = x-process continuous-time system matrix
%   fy = y-process continuous-time system matrix
%   fz = z-process continuous-time system matrix
%   put the F matrix from each process into one big F matrix. It will be a
%   diagonal matrix with a lot of zeros off-diagonal

if nargin == 2;
    fz = [];
end

lx1 = size(fx,1);
lx2 = size(fx,2);
ly1 = size(fy,1);
ly2 = size(fy,2);
lz1 = size(fz,1);
lz2 = size(fz,2);
% lx = length(fx);
% ly = length(fy);
% lz = length(fz);
%     n_states = lx+ly+lz;

F = [fx             zeros(lx1,ly2) zeros(lx1,lz2);
     zeros(ly1,lx2) fy             zeros(ly1,lz2);
     zeros(lz1,lx2) zeros(lz1,ly2) fz          ];

% if nargin == 2
%     lx = length(fx);
%     ly = length(fy);
% 
%     F = [fx           zeros(lx,ly);
%          zeros(ly,lx) fy          ];
end

