function [value,gradRho,gradY] = tracking(~,rho,y,cfg,parameters)

% Inputs
if isvector(rho)
    rho = rho(:);
end

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
value = weight.*cfg.dx.*(squaredDistance.'*rho);

% Gradients
gradRho = repmat(weight.*cfg.dx.*squaredDistance,1,size(rho,2));
gradY = zeros(size(y));

end
