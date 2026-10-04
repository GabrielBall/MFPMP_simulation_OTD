
meanFieldResults = run_simulation();

options.numParticles = 500;
options.numRealizations = 1;
options.comparisonEvery = 10;
options.initialSampling = 'quantile';
options.numComparisonCells = 50;
options.maxOverlayParticles = 500;

validation = n_system_validation(meanFieldResults,options);