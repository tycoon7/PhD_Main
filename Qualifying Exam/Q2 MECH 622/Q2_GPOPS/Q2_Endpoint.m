%-------------------------------------------%
% BEGIN: function Q4_Endpoint.m            % 
%-------------------------------------------%
function output = Q2_Endpoint(input)
% The possible outputs are "objective" and "eventgroup" I don't know what
% "eventgroup" is, but I'm pretty sure we are not using it.
% sprintf('Entering Endpoint')

T = input.phase.finaltime;

output.objective = T;

% sprintf('Leaving Endpoint')
end
%-------------------------------------------%
% END: function Q4_Endpoint.m              %
%-------------------------------------------%