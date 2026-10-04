function weight = gaussian(r,parameters)

% Parameters
b = parameters.bandwidth;
d = parameters.dimension;

% Normalisation
if parameters.normalised
    scale = 1/(2*pi*b^2)^(d/2);
else
    scale = 1;
end

% Kernel
weight = scale.*exp(-r.^2./(2*b^2));

end
