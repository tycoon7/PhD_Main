function [ m sigma ] = calcEnsembleStats( data )
%CALCENSEMBLESTATS returns the ensemble ave and std calculated from samples
%   data = input data is a matrix of realization data. The row dimension is
%          time and the column dimension is realizations
%   m = vector of the ensemble mean value at each time step
%   std = vector of the ensemble std value at each time step
%
%   To pass in one layer of a 3D matrix, use the squeeze command to isolate
%   the layer, then transpose as necessary

m = mean(data,2);     % calculated mean value of samples at each time
sigma = std(data,0,2);  % calculated std of samples at each time step

end

