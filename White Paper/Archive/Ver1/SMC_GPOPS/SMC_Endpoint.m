%-------------------------------------------%
% BEGIN: function SMC_Endpoint.m            % 
%-------------------------------------------%
function output = SMC_Endpoint(input)
% The possible outputs are "objective" and "eventgroup" I don't know what
% "eventgroup" is, but I'm pretty sure we are not using it.
% sprintf('Entering Endpoint')

S = input.phase.integral(1); % This is the integrated spot size

output.objective = S;

% sprintf('Leaving Endpoint')
end
%-------------------------------------------%
% END: function SMC_Endpoint.m              %
%-------------------------------------------%