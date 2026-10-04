function results = run_simulation()

% Setup
cfg = config();
rng(cfg.randomSeed);
rho0 = initialise_density(cfg);
y0 = cfg.initialLeader;

if cfg.sigmaY > 0
    cfg.leaderBrownianIncrements = sqrt(cfg.dt).*randn(1,cfg.numTimeSteps);
else
    cfg.leaderBrownianIncrements = zeros(1,cfg.numTimeSteps);
end

u = initialise_control(cfg);
cfg.interactionMatrix = build_interaction_matrix(cfg);

fprintf('Mean-field solver: boundary=%s, cells=%d, time steps=%d, sigmaX=%.3g\n',cfg.domain.boundary,cfg.numCells,cfg.numTimeSteps,cfg.sigmaX);

% Simulation
[J,Y,Rho] = simulate_forward(u,rho0,y0,cfg);
plot_solution(cfg,Y,Rho,u);

% Results
results.cfg = cfg;
results.objective = J;
results.control = u;
results.leader = Y;
results.density = Rho;

end

function u = initialise_control(cfg)

M = cfg.numTimeSteps;

switch lower(char(cfg.prescribedControl.type))
    case 'zero'
        u = zeros(1,M);

    case 'constant'
        u = cfg.prescribedControl.value.*ones(1,M);

    case 'function'
        if isempty(cfg.prescribedControl.functionHandle)
            error('No prescribed control function supplied.');
        end
        u = cfg.prescribedControl.functionHandle(cfg.time(1:end-1));
        u = reshape(u,1,[]);

    otherwise
        error('Unknown prescribed control type: %s.',cfg.prescribedControl.type);
end

if ~isequal(size(u),[1,M])
    error('The prescribed control has the wrong size.');
end

u = min(max(u,cfg.controlLower),cfg.controlUpper);

end
