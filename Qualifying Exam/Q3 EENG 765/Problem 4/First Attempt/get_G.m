function [ G ] = get_G( gx, gy, gz, varargin )
%GET_G Summary builds the big G matrix
%   Gx = x-process continuous-time noise input matrix
%   Gy = y-process continuous-time noise input matrix
%   Gz = z-process continuous-time noise input matrix
%   put the G matrix from each process into one big G matrix. It will be a
%   diagonal matrix with a lot of zeros off-diagonal

if nargin == 2;
    gz = [];
end

lx1 = size(gx,1);
lx2 = size(gx,2);
ly1 = size(gy,1);
ly2 = size(gy,2);
lz1 = size(gz,1);
lz2 = size(gz,2);

G = [gx             zeros(lx1,ly2) zeros(lx1,lz2);
     zeros(ly1,lx2) gy             zeros(ly1,lz2);
     zeros(lz1,lx2) zeros(lz1,ly2) gz          ];

end

