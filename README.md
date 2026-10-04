# Mean-Field Optimal Control for the Hegselmann--Krause Model

MATLAB implementation of a one-dimensional mean-field optimal-control framework for the Hegselmann--Krause opinion-dynamics model with a controlled leader.

## Repository structure

```text
MFPMP_simulation_OTD/
├── mean_field_PMP_simulation.prj
├── toolbox/
│   ├── config.m
│   ├── run_simulation.m
│   ├── simulate_forward.m
│   ├── solve_adjoints.m
│   ├── optimise_control.m
│   ├── tanh_continuation.m
│   ├── n_system_validation.m
│   ├── initialise_density.m
│   ├── forward_time_discretisation.m
│   ├── backward_time_discretisation.m
│   ├── plot_solution.m
│   ├── plot_optimisation_results.m
│   ├── +influence/
│   │   ├── +follower/
│   │   └── +leader/
│   ├── +loss/
│   │   ├── +control/
│   │   ├── +running/
│   │   └── +terminal/
│   ├── private/
│   └── examples/
└── resources/
```

## Getting started

The repository is self-contained MATLAB code and does not include third-party dependencies. Open `mean_field_PMP_simulation.prj` in MATLAB, or add the main toolbox and examples directory to the MATLAB path:

```matlab
addpath('toolbox')
addpath(fullfile('toolbox','examples'))
```

All main model and numerical parameters are defined in `toolbox/config.m`.

### Forward simulation

Run the default controlled mean-field simulation with

```matlab
results = run_simulation();
```

This initialises the density and prescribed control from `config.m`, solves the forward system, and plots the density evolution, final state, leader trajectory and control.

### Optimal control

A minimal optimisation is

```matlab
cfg = config();

rho0 = initialise_density(cfg);
y0 = cfg.initialLeader;
u0 = zeros(1,cfg.numTimeSteps);

results = optimise_control(u0,rho0,y0,cfg);

plot_optimisation_results( ...
    cfg, ...
    results.leader, ...
    results.density, ...
    results.control, ...
    results.costHistory, ...
    results.residualHistory);
```

The example `toolbox/examples/optimisation_script.m` performs the same basic workflow.

The returned structure contains the optimal control, density and leader trajectories, follower and leader adjoints, Hamiltonian gradient, active set, PMP residual, optimisation history and solver status.

### Kernel continuation

For a decreasing sequence of smooth transition widths,

```matlab
cfg = config();

rho0 = initialise_density(cfg);
y0 = cfg.initialLeader;
u0 = zeros(1,cfg.numTimeSteps);

lambdaValues = logspace(-1,-4,12);

continuation = tanh_continuation( ...
    lambdaValues, ...
    0, ...
    u0, ...
    rho0, ...
    y0, ...
    cfg);
```

Each stage is warm-started from the previous control.

### Particle-system validation

`n_system_validation.m` compares a mean-field solution against finite-\(N\) particle simulations. For example,

```matlab
meanFieldResults = run_simulation();

options.numParticles = 500;
options.numRealizations = 1;
options.comparisonEvery = 10;
options.initialSampling = 'quantile';
options.numComparisonCells = 50;
options.maxOverlayParticles = 500;

validation = n_system_validation(meanFieldResults,options);
```

The validation routine reports density errors and compares the empirical particle distribution with the mean-field solution over time. See `toolbox/examples/val_script.m` for a complete example.

A separate leader-free Hegselmann--Krause particle simulation is provided by `toolbox/examples/simulate_basic_hk.m`, with `plot_hk_mean_field.m` available for direct particle/mean-field comparison using the exact HK interaction rule.

## Included kernels and costs

Follower and leader influence functions are stored separately in `toolbox/+influence/`. The repository currently includes:

- Hegselmann--Krause indicator kernel
- tanh bounded-confidence approximation
- Gaussian kernel
- raised-cosine kernel
- compactly supported \(C^\infty\) bump kernel

For adjoint-based optimisation, a smooth leader kernel should be used. The exact HK indicator is primarily useful for forward simulation and comparison with the sharp bounded-confidence model.

The available loss functions are:

- quadratic control penalty;
- mean-square tracking running cost;
- binary-voter running cost;
- mean-square target terminal cost;
- binary-voter terminal cost.
