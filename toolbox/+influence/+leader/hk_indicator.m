function [weight,dWeightDr] = hk_indicator(r,parameters)

% Parameters
R = parameters.bandwidth;
alpha = parameters.alpha;

% Kernel
weight = alpha.*double(r <= R);
dWeightDr = zeros(size(r));

end
