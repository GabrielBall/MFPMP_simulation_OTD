function continuation = tanh_continuation(lambdaValues,finalLambda,u0,rho0,y0,cfg)

%% Input Gathering
lambdaValues = double(lambdaValues(:));
numberOfStages = numel(lambdaValues);
M = cfg.numTimeSteps;
Nx = cfg.numCells;

rho0 = rho0(:);
%% Input Validation

if numberOfStages < 2
    error('At least two lambda values are required.');
end

if any(~isfinite(lambdaValues)) || any(lambdaValues <= 0)

    error('All lambda values must be positive and finite.');
end

if any(diff(lambdaValues) >= 0)
    error('lambdaValues must be strictly decreasing.');
end

if ~isequal(size(u0),[1,M])
    error('u0 must have size 1-by-cfg.numTimeSteps.');
end

if ~isscalar(y0)
    error('y0 must be scalar.');
end

if finalLambda ~= 0
    error('finalLambda must be zero. The final stage is an exact-HK forward simulation, not another smooth optimisation.');
end

%% Allocation
controlHistory =zeros(numberOfStages,M);
leaderHistory = zeros(numberOfStages,M+1);
objectiveHistory = zeros(numberOfStages,1);
controlDifference = nan(numberOfStages,1);
pmpResidualHistory = nan(numberOfStages,1);
leaderDifference = nan(numberOfStages,1);
objectiveDifference = nan(numberOfStages,1);
exitflagHistory = nan(numberOfStages,1);
statusHistory = cell(numberOfStages,1);


warmStartControl = u0;

fprintf('\n');
fprintf('HK kernel continuation\n');
fprintf('Number of stages: %d\n',numberOfStages);

if ~isfield(cfg,'interactionMatrix') || isempty(cfg.interactionMatrix)
    cfg.interactionMatrix = build_interaction_matrix(cfg);
end


%% Continuation loop
for stage = 1:numberOfStages

    lambda = lambdaValues(stage);
    cfgStage = cfg;
    cfgStage.leaderKernelParameters.transitionWidth = lambda;

    fprintf('\n');
    fprintf('Continuation stage %d/%d: lambda = %.12e\n',stage,numberOfStages,lambda);

    results = optimise_control(warmStartControl,rho0,y0,cfgStage);

    controlHistory(stage,:) =results.control;
    leaderHistory(stage,:) =results.leader;
    objectiveHistory(stage) =results.objective;
    pmpResidualHistory(stage) =results.pmpResidual;
    statusHistory{stage} =results.status;

    if isfield(results,'exitflag')
        exitflagHistory(stage) = results.exitflag;
    end

    if stage > 1

        controlIncrement = controlHistory(stage,:)-controlHistory(stage-1,:);
        leaderIncrement = leaderHistory(stage,:) - leaderHistory(stage-1,:);

        controlDifference(stage) = sqrt(cfgStage.dt) .*norm(controlIncrement,2);
        leaderDifference(stage) = max(abs(leaderIncrement));
        objectiveDifference(stage) = abs(objectiveHistory(stage) - objectiveHistory(stage-1));

        fprintf('U difference: %.12e\n',controlDifference(stage));
        fprintf('Y difference: %.12e\n',leaderDifference(stage));
        fprintf('J difference: %.12e\n',objectiveDifference(stage));
    end

    fprintf('Objective: %.12e\n',objectiveHistory(stage));
    fprintf('PMP residual: %.12e\n',pmpResidualHistory(stage));
    fprintf('Status: %s\n',statusHistory{stage});
    
    warmStartControl =results.control;
end

cfgFinal = cfg;

cfgFinal.leaderKernel = @influence.leader.hk_indicator;

fprintf('\n');
fprintf('Final exact-HK forward simulation\n');

[JFinal,YFinal,RhoFinal] = simulate_forward(warmStartControl,rho0,y0,cfgFinal);

finalResults.cfg = cfgFinal;
finalResults.status = 'exact_hk_forward_simulation';
finalResults.exitflag = NaN;
finalResults.objective = JFinal;
finalResults.control = warmStartControl;
finalResults.leader = YFinal;
finalResults.density = RhoFinal;
finalResults.costHistory = JFinal;
finalResults.residualHistory = [];
finalResults.stepHistory = [];

%% Data Packaging

continuation.lambda = lambdaValues;
continuation.control = controlHistory;
continuation.leader = leaderHistory;
continuation.objective = objectiveHistory;
continuation.controlDifference = controlDifference;
continuation.leaderDifference = leaderDifference;
continuation.objectiveDifference = objectiveDifference;
continuation.pmpResidual = pmpResidualHistory;
continuation.exitflag = exitflagHistory;
continuation.status = statusHistory;
continuation.finalResults = finalResults;
continuation.finalLambda = finalLambda;

continuation.summary = table( ...
        lambdaValues, ...
        objectiveHistory, ...
        controlDifference, ...
        leaderDifference, ...
        objectiveDifference, ...
        pmpResidualHistory, ...
        exitflagHistory, ...
        statusHistory, ...
        'VariableNames',{ ...
            'lambda', ...
            'objective', ...
            'controlDifference', ...
            'leaderDifference', ...
            'objectiveDifference', ...
            'pmpResidual', ...
            'exitflag', ...
            'status'});

%% Plot convergence against lambda

plot_hk_continuation(lambdaValues,controlHistory,leaderHistory,objectiveHistory,cfg);

%% Plot the final exact-HK system

fprintf('\n');
fprintf('Plotting the final exact-HK system\n');

plot_optimisation_results(finalResults.cfg,finalResults.leader,finalResults.density,finalResults.control,finalResults.costHistory,finalResults.residualHistory);

end

function plot_hk_continuation(lambdaValues,controlHistory,leaderHistory,objectiveHistory,cfg)

legendEntries = arrayfun(@(lambda) sprintf('\\lambda = %.3g',lambda),lambdaValues(:),'UniformOutput',false);

figure('Name','Continuation controls');
stairs(cfg.time(1:end-1),controlHistory.','LineWidth',1.5);
xlabel('t');
ylabel('u(t)');
legend(legendEntries,'Location','best');
grid on;

figure('Name','Continuation leader trajectories');
plot(cfg.time,leaderHistory.','LineWidth',1.5);
hold on;
yline(cfg.target,'--','LineWidth',1.5);
xlabel('t');
ylabel('y(t)');
legend([legendEntries;{'Target'}],'Location','best');
grid on;

figure('Name','Continuation objective');
semilogx(lambdaValues,objectiveHistory,'o-','LineWidth',1.5);
set(gca,'XDir','reverse');
xlabel('\lambda');
ylabel('J');
grid on;

end
