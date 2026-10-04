function [weight,dWeightDr] = tanh_approx(r,parameters)

% Parameters
alpha = parameters.alpha;
R = parameters.bandwidth;
delta = parameters.transitionWidth;

% Kernel
argument = (R.^2-r.^2)./(2.*R.*delta);
tanhArgument = tanh(argument);
weight = 0.5.*alpha.*(1+tanhArgument);
dWeightDr = -alpha.*r./(2.*R.*delta).*(1-tanhArgument.^2);

end
