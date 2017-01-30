%-------------------------------------------%
% BEGIN: function Q4_Endpoint.m            % 
%-------------------------------------------%
function output = Q4_Endpoint(input)
% The possible outputs are "objective" and "eventgroup" I don't know what
% "eventgroup" is, but I'm pretty sure we are not using it.
% sprintf('Entering Endpoint')

E = input.phase.integral; % This is the expended energy

output.objective = E;

% sprintf('Leaving Endpoint')
end
%-------------------------------------------%
% END: function Q4_Endpoint.m              %
%-------------------------------------------%