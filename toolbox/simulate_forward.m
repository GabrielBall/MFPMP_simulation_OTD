function [J,y,rho] = simulate_forward(u,rho0,y0,cfg)

% Inputs
Nt = cfg.numTimeSteps;
Nx = cfg.numCells;
dt = cfg.dt;
rho0 = rho0(:);
boundaryCondition = validatestring(lower(char(cfg.domain.boundary)),{'periodic','zero_flux','open'});

if isfield(cfg,'interactionMatrix') && ~isempty(cfg.interactionMatrix)
    interactionMatrix = cfg.interactionMatrix;
else
    interactionMatrix = build_interaction_matrix(cfg);
end

if isfield(cfg,'leaderBrownianIncrements')
    leaderBrownianIncrements = cfg.leaderBrownianIncrements;
else
    leaderBrownianIncrements = zeros(1,Nt);
end

% Validation
if ~isequal(size(u),[1,Nt])
    error('u must have size 1-by-numTimeSteps.');
end
if numel(rho0) ~= Nx
    error('rho0 must contain one value per cell.');
end
if ~isequal(size(interactionMatrix),[Nx+1,Nx])
    error('The interaction matrix has the wrong size.');
end

% Allocation
y = zeros(1,Nt+1);
rho = zeros(Nx,Nt+1);
runningLoss = zeros(1,Nt);
y(1) = apply_leader_domain(y0,cfg.domain,boundaryCondition);
rho(:,1) = rho0;
A = @(rho,y) chang_cooper_transport_operator(rho,y,cfg,interactionMatrix,boundaryCondition);

% Forward loop
for k = 1:Nt
    rhoPrev = rho(:,k);
    yPrev = y(k);

    leaderNoise = cfg.sigmaY.*leaderBrownianIncrements(k);
    yNext = yPrev+dt.*u(k)+leaderNoise;
    y(k+1) = apply_leader_domain(yNext,cfg.domain,boundaryCondition);

    rho(:,k+1) = forward_time_discretisation(rhoPrev,dt,cfg.timescheme,A,[],{yPrev},{y(k+1)});

    if min(rho(:,k+1)) < -cfg.numerics.positivityTolerance
        error('Density positivity failure at step %d: minimum %.3e.',k,min(rho(:,k+1)));
    end

    runningStateCost = cfg.runningCost(cfg.time(k),rhoPrev,yPrev,cfg,cfg.runningCostParameters);
    controlCost = cfg.controlPenalty(cfg.time(k),u(k),cfg,cfg.controlPenaltyParameters);
    runningLoss(k) = runningStateCost+controlCost;
end

% Objective
terminalCost = cfg.terminalCost(rho(:,end),y(end),cfg,cfg.terminalCostParameters);
J = dt.*sum(runningLoss)+terminalCost;

end

function state = apply_leader_domain(state,domain,boundaryCondition)

if strcmp(boundaryCondition,'periodic')
    domainWidth = domain.upper-domain.lower;
    state = domain.lower+mod(state-domain.lower,domainWidth);
end

end