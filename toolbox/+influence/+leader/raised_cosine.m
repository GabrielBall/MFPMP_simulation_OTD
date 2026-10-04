function [weight,dWeightDr] = raised_cosine(r,parameters)

% Parameters
alpha = parameters.alpha;
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
dWeightDr = zeros(size(r));
weight(core) = alpha;

if any(taper(:))
    theta = pi.*(r(taper)-coreRadius)./delta;
    weight(taper) = 0.5.*alpha.*(1+cos(theta));
    dWeightDr(taper) = -alpha.*pi./(2.*delta).*sin(theta);
end

end
