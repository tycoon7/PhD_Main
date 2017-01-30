%-------------------------------------------%
% BEGIN: function Q5P2_Endpoint.m            % 
%-------------------------------------------%
function output = Q5P2_Endpoint(input)
% The possible outputs are "objective" and "eventgroup" I don't know what
% "eventgroup" is, but I'm pretty sure we are not using it.
% sprintf('Entering Endpoint')

R = input.parameter;        % Radius of mirror
J = input.phase.integral;   % This is the integrated radius

% Penalize shorter focal length and focal point evolution ball radius
output.objective = (1-R) + J;   

% global radius
% output.objective = radius;

% sprintf('Leaving Endpoint')
end
%-------------------------------------------%
% END: function Q5P2_Endpoint.m              %
%-------------------------------------------%