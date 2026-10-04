function results = optimise_control(u0,rho0,y0,cfg)

%% Input gathering
M = cfg.numTimeSteps;
Nx = cfg.numCells;
dt = cfg.dt;
rho0 = rho0(:);

optimisation = cfg.optimisation;
maxIterations = optimisation.maxIterations;
pmpTolerance = optimisation.gradientTolerance;
maxTrials = optimisation.maxLineSearches;

sqhEpsilon = optimisation.sqhEpsilon;
sqhIncrease = optimisation.sqhIncrease;
sqhDecrease = optimisation.sqhDecrease;
sqhEta = optimisation.sqhEta;
boundTolerance = optimisation.sqhBoundTolerance;

lambda = cfg.controlPenaltyParameters.lambda;

%% Input validation
if ~isequal(size(u0),[1,M])
    error('u0 must have size 1-by-cfg.numTimeSteps.');
end
if numel(rho0) ~= Nx
    error('rho0 must contain cfg.numCells entries.');
end
if ~isscalar(y0)
    error('y0 must be scalar.');
end
if strcmpi(cfg.domain.boundary,'periodic')
    error('optimise_control requires non-periodic leader dynamics.');
end
if get_option(optimisation,'enforceLeaderStateBounds',false)
    error('optimise_control does not support an absolute leader-state constraint.');
end
if sqhEpsilon <= 0 || sqhIncrease <= 1 || sqhDecrease <= 0 || sqhDecrease >= 1 || sqhEta <= 0
    error('The SQH parameters are invalid.');
end

controlLower = expand_bound(cfg.controlLower,[1,M]);
controlUpper = expand_bound(cfg.controlUpper,[1,M]);
if any(controlLower(:) > controlUpper(:))
    error('A lower control bound exceeds an upper control bound.');
end

if ~isfield(cfg,'interactionMatrix') || isempty(cfg.interactionMatrix)
    cfg.interactionMatrix = build_interaction_matrix(cfg);
end

u = min(max(u0,controlLower),controlUpper);

%% Initial solve
[J,y,rho] = simulate_forward(u,rho0,y0,cfg);
[P,Q] = solve_adjoints(y,rho,cfg);
G = hamiltonian_gradient(u,Q,cfg);
[restrictedG,active] = restricted_hamiltonian_gradient(u,G,controlLower,controlUpper,boundTolerance);
pmpResidual = norm(restrictedG(:))./sqrt(M);

costHistory = J;
residualHistory = pmpResidual;
epsilonHistory = sqhEpsilon;
updateHistory = [];
forwardSolves = 1;
adjointSolves = 1;
acceptedIterations = 0;
status = 'maximum_iterations';

%% Main loop
for iteration = 1:maxIterations
    fprintf('Iteration %d: J = %.10e, PMP residual = %.3e, epsilon = %.3e\n',iteration,J,pmpResidual,sqhEpsilon);

    if pmpResidual <= pmpTolerance
        status = 'converged';
        break;
    end

    accepted = false;

    for trial = 1:maxTrials
        uTrial = (2.*sqhEpsilon.*u-Q(2:end))./(lambda+2.*sqhEpsilon);
        uTrial = min(max(uTrial,controlLower),controlUpper);

        tau = dt.*sum((uTrial-u).^2);

        if tau == 0
            status = 'zero_sqh_step';
            break;
        end

        [JTrial,yTrial,rhoTrial] = simulate_forward(uTrial,rho0,y0,cfg);
        forwardSolves = forwardSolves+1;

        if JTrial-J <= -sqhEta.*tau
            accepted = true;
            break;
        end

        sqhEpsilon = sqhIncrease.*sqhEpsilon;
    end

    if ~accepted
        if ~strcmp(status,'zero_sqh_step')
            status = 'sqh_failed';
        end
        break;
    end

    u = uTrial;
    J = JTrial;
    y = yTrial;
    rho = rhoTrial;

    [P,Q] = solve_adjoints(y,rho,cfg);
    adjointSolves = adjointSolves+1;

    G = hamiltonian_gradient(u,Q,cfg);
    [restrictedG,active] = restricted_hamiltonian_gradient(u,G,controlLower,controlUpper,boundTolerance);
    pmpResidual = norm(restrictedG(:))./sqrt(M);

    acceptedIterations = acceptedIterations+1;
    costHistory(end+1,1) = J;
    residualHistory(end+1,1) = pmpResidual;
    updateHistory(end+1,1) = tau;

    sqhEpsilon = sqhDecrease.*sqhEpsilon;
    epsilonHistory(end+1,1) = sqhEpsilon;
end

if pmpResidual <= pmpTolerance
    status = 'converged';
end

results.cfg = cfg;
results.status = status;
results.objective = J;
results.control = u;
results.leader = y;
results.density = rho;
results.densityAdjoint = P;
results.leaderAdjoint = Q;
results.gradient = G;
results.restrictedGradient = restrictedG;
results.activeSet = active;
results.pmpResidual = pmpResidual;
results.costHistory = costHistory;
results.residualHistory = residualHistory;
results.epsilonHistory = epsilonHistory;
results.updateHistory = updateHistory;
results.optimisationOutput.algorithm = 'sequential quadratic Hamiltonian';
results.optimisationOutput.iterations = acceptedIterations;
results.optimisationOutput.forwardSolves = forwardSolves;
results.optimisationOutput.adjointSolves = adjointSolves;

end

function G = hamiltonian_gradient(u,Q,cfg)
M = cfg.numTimeSteps;
if ~isequal(size(Q),[1,M+1])
    error('Q has the wrong size.');
end

[~,controlGradient] = cfg.controlPenalty(cfg.time(1:end-1),u,cfg,cfg.controlPenaltyParameters);
controlGradient = reshape(controlGradient,1,M);
G = Q(2:end)+controlGradient;

if any(~isfinite(G))
    error('The Hamiltonian gradient is non-finite.');
end
end

function [restrictedG,active] = restricted_hamiltonian_gradient(u,G,controlLower,controlUpper,boundTolerance)
atLower = u <= controlLower+boundTolerance;
atUpper = u >= controlUpper-boundTolerance;
active = (atLower & G >= 0) | (atUpper & G <= 0);
restrictedG = G;
restrictedG(active) = 0;
end

function expanded = expand_bound(value,controlSize)
if isscalar(value)
    expanded = repmat(value,controlSize);
elseif isequal(size(value),controlSize)
    expanded = value;
elseif isvector(value) && numel(value) == controlSize(2)
    expanded = reshape(value,controlSize);
else
    error('A control bound has the wrong size.');
end
end

function value = get_option(options,fieldName,defaultValue)
if isfield(options,fieldName)
    value = options.(fieldName);
else
    value = defaultValue;
end
end
