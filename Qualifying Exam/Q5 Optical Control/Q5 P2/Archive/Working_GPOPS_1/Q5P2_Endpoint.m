%-------------------------------------------%
% BEGIN: function Q5P2_Endpoint.m            % 
%-------------------------------------------%
function output = Q5P2_Endpoint(input)
% The possible outputs are "objective" and "eventgroup" I don't know what
% "eventgroup" is, but I'm pretty sure we are not using it.
% sprintf('Entering Endpoint')

J = input.phase.integral;           % This is the cost function

output.objective = J;

% sprintf('Leaving Endpoint')
end
%-------------------------------------------%
% END: function Q5P2_Endpoint.m              %
%-------------------------------------------%