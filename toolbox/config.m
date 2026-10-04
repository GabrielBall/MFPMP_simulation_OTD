function cfg = config()
%% General Settings
cfg.randomSeed   = 42;
cfg.sigmaX = 0.05;
cfg.sigmaY = 0.0;

%% Time Discretisation
cfg.dimension    = 1;
cfg.finalTime    = 100.0;
cfg.numTimeSteps = 1000;
cfg.dt           = cfg.finalTime/cfg.numTimeSteps;
cfg.time         = linspace(0,cfg.finalTime,cfg.numTimeSteps+1);
cfg.timescheme   = 'semi_implicit';

%% Domain Setup
boundaryCondition = 'zero_flux';
boundaryCondition = validatestring(lower(char(boundaryCondition)),{'periodic','zero_flux','open'});

cfg.domain.type     = 'interval';
cfg.domain.lower    = -1.0;
cfg.domain.upper    = 1.0;
cfg.domain.boundary = boundaryCondition;

%% Domain Discretisation
cfg.numCells = 100;
cfg.dx =(cfg.domain.upper-cfg.domain.lower)/cfg.numCells;

cfg.cellEdges = linspace(cfg.domain.lower, cfg.domain.upper, cfg.numCells+1).';

cfg.cellCentres = 0.5 .* (cfg.cellEdges(1:end-1) + cfg.cellEdges(2:end));

%% Initial Conditions
cfg.initialDensity.type = 'uniform';
cfg.initialDensity.lower = -1.0;
cfg.initialDensity.upper = 1.0;
cfg.initialDensity.mean = 0.0;
cfg.initialDensity.standardDeviation = 0.35;
cfg.initialDensity.values = [];
cfg.initialDensity.functionHandle = [];

%% Follower Influence
cfg.followerKernel = @influence.follower.tanh_approx;

cfg.followerKernelParameters.bandwidth = 0.3;
cfg.followerKernelParameters.dimension = 1;
cfg.followerKernelParameters.transitionWidth = 0.001;
cfg.followerKernelParameters.normalised = false;

%% Leader Influence and Target
cfg.leaderKernel = @influence.leader.tanh_approx;

cfg.leaderKernelParameters.alpha = 0.05;
cfg.leaderKernelParameters.bandwidth = 0.3;
cfg.leaderKernelParameters.dimension = 1;
cfg.leaderKernelParameters.normalised = false;
cfg.leaderKernelParameters.transitionWidth = 0.01;

cfg.leaderFieldScale = 1.0;

cfg.initialLeader = -0.0;
cfg.target        = -0.5;

cfg.controlLower = -0.1;
cfg.controlUpper = 0.1;

cfg.prescribedControl.type = 'zero';
cfg.prescribedControl.value = 0.0;
cfg.prescribedControl.functionHandle = [];
%% Cost Parameters
cfg.controlPenalty = @loss.control.quadratic;

cfg.controlPenaltyParameters.lambda = 0.01;

cfg.runningCost = @loss.running.tracking;
cfg.runningCostParameters.weight = 0.0;
cfg.runningCostParameters.leaderWeight = 1.0;
cfg.runningCostParameters.target = cfg.target;

cfg.terminalCost = @loss.terminal.mean_square_target;

cfg.terminalCostParameters.weight = 1.0;
cfg.terminalCostParameters.target = cfg.target;

cfg.terminalCostParameters.threshold = 0;
cfg.terminalCostParameters.transitionWidth = 0.1;

%% Optimisation Parameters
cfg.optimisation.maxIterations = 1000;
cfg.optimisation.gradientTolerance = 1e-4;
cfg.optimisation.sqhMaxTrials = 30;

cfg.optimisation.sqhBoundTolerance = 1e-10;

cfg.optimisation.sqhEpsilon = 100;
cfg.optimisation.sqhIncrease = 25;
cfg.optimisation.sqhDecrease = 0.1;
cfg.optimisation.sqhEta = 1e-8;
cfg.optimisation.sqhKappa = 0;
cfg.optimisation.sqhControlWeight = 1;
cfg.optimisation.sqhLeaderWeight = 1;

%% Misc
cfg.numerics.positivityTolerance = 1e-11;
cfg.numerics.interactionStorage = 'auto';

end
