function weight = raised_cosine(r,parameters)

% Parameters
R = parameters.bandwidth;
delta = 0.2*R;

if isfield(parameters,'transitionWidth')
    delta = parameters.transitionWidth;
end

% Regions
coreRadius = R-delta;
core = r <= coreRadius;
taper = r > coreRadius & r < R;

% Kernel
weight = zeros(size(r));
weight(core) = 1;

if any(taper(:))
    theta = pi.*(r(taper)-coreRadius)./delta;
    weight(taper) = 0.5.*(1+cos(theta));
end

end
