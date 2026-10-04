function rho0 = initialise_density(cfg)

% Inputs
edges = cfg.cellEdges(:);
centres = cfg.cellCentres(:);
dx = cfg.dx;
type = lower(char(cfg.initialDensity.type));

% Initialisation
switch type
    case 'uniform'
        if cfg.initialDensity.upper <= cfg.initialDensity.lower
            error('The uniform density interval must have positive width.');
        end

        lowerMembership = max(edges(1:end-1),cfg.initialDensity.lower);
        upperMembership = min(edges(2:end),cfg.initialDensity.upper);
        membership = max(upperMembership-lowerMembership,0);
        rho0 = membership./((cfg.initialDensity.upper-cfg.initialDensity.lower).*dx);

    case {'normal','gaussian'}
        scaledEdges = (edges-cfg.initialDensity.mean)./(sqrt(2).*cfg.initialDensity.standardDeviation);
        mass = diff(0.5.*(1+erf(scaledEdges)));
        totalMass = sum(mass);

        if totalMass <= 0
            error('The Gaussian has zero numerical mass on the domain.');
        end

        rho0 = mass./(totalMass.*dx);

    case 'specified'
        rho0 = cfg.initialDensity.values(:);
        if numel(rho0) ~= cfg.numCells
            error('initialDensity.values must contain one value per spatial cell.');
        end

    case 'function'
        f = cfg.initialDensity.functionHandle;
        if isempty(f)
            error('initialDensity.functionHandle must be supplied for a function initial density.');
        end
        rho0 = f(centres);
        rho0 = rho0(:);

    otherwise
        error('Unknown initial density type: %s.',cfg.initialDensity.type);
end


% Checks
if any(~isfinite(rho0)) || any(rho0 < 0)
    error('The initial density must be finite and non-negative.');
end

mass = dx.*sum(rho0);
if mass <= 0
    error('The initial density has zero mass.');
end
rho0 = rho0./mass;

end
