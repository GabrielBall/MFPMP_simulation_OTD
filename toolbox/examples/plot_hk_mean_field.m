function comparison = plot_hk_mean_field(numParticles,sigmaX)
%PLOT_HK_MEAN_FIELD Compare the leader-free HK particle system and mean field.
%
%   comparison = plot_hk_mean_field(numParticles,sigmaX)
%
% The particle system and mean-field equation use the same exact HK kernel,
% initial distribution, noise strength, time grid and boundary condition.

cfg = config();

if nargin < 1 || isempty(numParticles)
    numParticles = 500;
end

if nargin < 2 || isempty(sigmaX)
    sigmaX = cfg.sigmaX;
end

particleResults = simulate_basic_hk(numParticles,sigmaX);

cfg.sigmaX = sigmaX;
cfg.followerKernel = @hk_indicator;
cfg.leaderKernel = @zero_kernel;
cfg.leaderKernelParameters = struct();
cfg.leaderFieldScale = 0;
cfg.sigmaY = 0;
cfg.runningCost = @zero_running_cost;
cfg.runningCostParameters = struct();
cfg.controlPenalty = @zero_control_cost;
cfg.controlPenaltyParameters = struct();
cfg.terminalCost = @zero_terminal_cost;
cfg.terminalCostParameters = struct();

rho0 = initialise_density(cfg);
u = zeros(1,cfg.numTimeSteps);
[~,~,rho] = simulate_forward(u,rho0,0,cfg);

figure('Name','Leader-free Hegselmann-Krause: particles and mean field');
imagesc(cfg.time,cfg.cellCentres,rho);
axis xy tight;
hold on;

maxTrajectories = 200;
indices = unique(round(linspace(1,numParticles,min(numParticles,maxTrajectories))));
plot(cfg.time,particleResults.particles(indices,:).','k-','LineWidth',0.35);

xlabel('Time');
ylabel('Opinion');
title(sprintf('HK particle system and mean-field approximation, N = %d, \\sigma = %.3g',numParticles,sigmaX));
colorbar;
colormap(turbo);

comparison.cfg = cfg;
comparison.particles = particleResults.particles;
comparison.density = rho;
comparison.time = cfg.time;
comparison.sigmaX = sigmaX;
comparison.numParticles = numParticles;

end

function weight = hk_indicator(r,parameters)

R = parameters.bandwidth;
weight = double(r <= R);

end

function weight = zero_kernel(r,~)

weight = zeros(size(r));

end

function value = zero_running_cost(varargin)

value = 0;

end

function value = zero_control_cost(varargin)

value = 0;

end

function value = zero_terminal_cost(varargin)

value = 0;

end
