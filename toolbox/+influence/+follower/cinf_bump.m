function weight = cinf_bump(r,parameters)

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
    t = (r(taper)-coreRadius)./delta;
    weight(taper) = smooth_cutoff(t);
end

end

function h = smooth_cutoff(t)

% Stable cutoff
z = 1./(1-t)-1./t;
h = zeros(size(t));
left = z <= 0;
h(left) = 1./(1+exp(z(left)));
h(~left) = exp(-z(~left))./(1+exp(-z(~left)));

end
