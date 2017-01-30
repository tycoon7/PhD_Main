function [ T ] = titles( t )
%TITLES this function will allow plot titles to be set in a loop
%   have to use cellstr to make a string an array entry

Tarray = cellstr(['Centered & Ideal'; 'Off-Axis & Ideal'; 'Off-Axis & Rough']);
              
T = Tarray(t);
end