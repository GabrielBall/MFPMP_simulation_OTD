function [value,gradRho,gradY] = mean_square_target(rho,y,cfg,parameters)

% Inputs
rho = rho(:);
target = cfg.target;

if isfield(parameters,'target')
    target = parameters.target;
end

% Distance
displacement = cfg.cellCentres(:)-target;

if strcmpi(cfg.domain.boundary,'periodic')
    width = cfg.domain.upper-cfg.domain.lower;
    displacement = mod(displacement+0.5.*width,width)-0.5.*width;
end

% Cost
weight = parameters.weight;
squaredDistance = displacement.^2;
value = weight.*cfg.dx.*sum(squaredDistance.*rho);

% Gradients
gradRho = weight.*cfg.dx.*squaredDistance;
gradY = zeros(size(y));

end
