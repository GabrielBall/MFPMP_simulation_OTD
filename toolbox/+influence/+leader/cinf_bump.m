function [weight,dWeightDr] = cinf_bump(r,parameters)

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
    t = (r(taper)-coreRadius)./delta;
    [h,dhdt] = smooth_cutoff(t);
    weight(taper) = alpha.*h;
    dWeightDr(taper) = alpha.*dhdt./delta;
end

end

function [h,dhdt] = smooth_cutoff(t)

% Stable cutoff
z = 1./(1-t)-1./t;
h = zeros(size(t));
left = z <= 0;
h(left) = 1./(1+exp(z(left)));
h(~left) = exp(-z(~left))./(1+exp(-z(~left)));

% Derivative
dhdt = zeros(size(t));
active = abs(z) < 350;

if any(active(:))
    logisticDerivative = h(active).*(1-h(active));
    zPrime = 1./t(active).^2+1./(1-t(active)).^2;
    dhdt(active) = -logisticDerivative.*zPrime;
end

end
