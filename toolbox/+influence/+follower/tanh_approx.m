function weight = tanh_approx(r,parameters)

% Parameters
R = parameters.bandwidth;
delta = parameters.transitionWidth;

% Kernel
argument = (R.^2-r.^2)./(2.*R.*delta);
weight = 0.5.*(1+tanh(argument));

end
