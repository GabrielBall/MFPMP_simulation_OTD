cfg = config();

rho0 = initialise_density(cfg);
y0 = cfg.initialLeader;
u0 = zeros(1,cfg.numTimeSteps);

lambdaValues = cfg.leaderKernelParameters.bandwidth .* logspace(0,-10,10);

finalLambda = 0;

continuation = tanh_continuation( ...
        lambdaValues, ...
        finalLambda, ... 
        u0, ...
        rho0, ...
        y0, ...
        cfg);


disp(continuation.summary);