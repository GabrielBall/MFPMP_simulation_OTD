cfg = config();

rho0 = initialise_density(cfg);
y0 = cfg.initialLeader;
u0 = zeros(1,cfg.numTimeSteps);

results = optimise_control(u0,rho0,y0,cfg);

plot_optimisation_results(cfg, ...
    results.leader ...
    ,results.density ...
    ,results.control)
