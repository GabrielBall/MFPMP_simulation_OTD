function [weight,dWeightDr] = gaussian(r,parameters)

% Parameters
alpha = parameters.alpha;
b = parameters.bandwidth;
d = parameters.dimension;

% Normalisation
if parameters.normalised
    scale = 1/(2*pi*b^2)^(d/2);
else
    scale = 1;
end

% Kernel
weight = alpha.*scale.*exp(-r.^2./(2*b^2));
dWeightDr = -(r./b^2).*weight;

end
