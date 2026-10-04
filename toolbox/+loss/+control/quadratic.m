function [value,gradU] = quadratic(~,u,~,parameters)

% Penalty
lambda = parameters.lambda;
value = 0.5.*lambda.*sum(u.^2,1);
gradU = lambda.*u;

end
