function results = simulate_basic_hk(numParticles,sigmaX)
%SIMULATE_BASIC_HK Simulate the leader-free continuous-time HK particle model.
%
%   results = simulate_basic_hk(numParticles,sigmaX)
%
% Simulates
%   dX_i = -(1/N) sum_j phi(|X_i-X_j|)(X_i-X_j) dt + sigmaX dW_i
% using the exact HK bounded-confidence kernel phi(r)=1_{r<=R}.
%
% The time grid, domain, confidence radius and initial distribution are taken
% from config.m. Set sigmaX = 0 for the deterministic model.

cfg = config();

if nargin < 1 || isempty(numParticles)
    numParticles = 500;
end

if nargin < 2 || isempty(sigmaX)
    sigmaX = cfg.sigmaX;
end

rng(cfg.randomSeed);

N = numParticles;
Nt = cfg.numTimeSteps;
dt = cfg.dt;
R = cfg.followerKernelParameters.bandwidth;

rho0 = initialise_density(cfg);
X = zeros(N,Nt+1);
X(:,1) = sample_initial_particles(rho0,cfg,N);

for k = 1:Nt
    x = X(:,k);
    displacement = x-x.';

    if strcmpi(cfg.domain.boundary,'periodic')
        width = cfg.domain.upper-cfg.domain.lower;
        displacement = mod(displacement+0.5*width,width)-0.5*width;
    end

    weight = double(abs(displacement) <= R);
    velocity = -sum(weight.*displacement,2)/N;

    x = x+dt*velocity+sigmaX*sqrt(dt)*randn(N,1);

    switch lower(cfg.domain.boundary)
        case 'periodic'
            width = cfg.domain.upper-cfg.domain.lower;
            x = cfg.domain.lower+mod(x-cfg.domain.lower,width);
        case 'zero_flux'
            width = cfg.domain.upper-cfg.domain.lower;
            shifted = mod(x-cfg.domain.lower,2*width);
            x = cfg.domain.lower+min(shifted,2*width-shifted);
        otherwise
            error('simulate_basic_hk supports periodic and zero_flux boundaries.');
    end

    X(:,k+1) = x;
end

results.cfg = cfg;
results.time = cfg.time;
results.particles = X;
results.numParticles = N;
results.sigmaX = sigmaX;
results.confidenceRadius = R;

end

function X = sample_initial_particles(rho0,cfg,N)

cellMass = max(rho0(:),0)*cfg.dx;
cellMass = cellMass/sum(cellMass);
cdf = cumsum(cellMass);
cdf(end) = 1;

u = rand(N,1);
cellIndex = discretize(u,[0;cdf]);
lowerCdf = [0;cdf(1:end-1)];
withinCell = (u-lowerCdf(cellIndex))./cellMass(cellIndex);
X = cfg.cellEdges(cellIndex)+cfg.dx*withinCell;

end
