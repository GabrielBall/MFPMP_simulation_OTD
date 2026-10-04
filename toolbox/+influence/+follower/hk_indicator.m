function weight = hk_indicator(r,parameters)

% Parameters
R = parameters.bandwidth;
alpha = 1;

if isfield(parameters,'alpha')
    alpha = parameters.alpha;
end

% Kernel
weight = alpha.*double(r <= R);

end
