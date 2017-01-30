function [ T ] = titles1( t )
%TITLES this function will allow plot titles to be set in a loop
%   have to use cellstr to make a string an array entry

Tarray = cellstr([  'Force Input, x_1 Output'; ... 
                    'Force Input, x_2 output';      ]);
              
T = Tarray(t);
end

