function interactionMatrix = build_interaction_matrix(cfg)

% Geometry
interfaceLocation = cfg.cellEdges(:);
cellLocation = cfg.cellCentres(:);
displacement = interfaceLocation-cellLocation.';
boundaryCondition = lower(char(cfg.domain.boundary));

if strcmp(boundaryCondition,'periodic')
    displacement = periodic_displacement(displacement,cfg.domain);
elseif ~strcmp(boundaryCondition,'zero_flux') && ~strcmp(boundaryCondition,'open')
    error('Unsupported boundary condition: %s.',cfg.domain.boundary);
end

% Interaction
followerWeight = cfg.followerKernel(abs(displacement),cfg.followerKernelParameters);
interactionMatrix = -cfg.dx.*followerWeight.*displacement;

% Storage
switch lower(char(cfg.numerics.interactionStorage))
    case 'dense'
    case 'sparse'
        interactionMatrix = sparse(interactionMatrix);
    case 'auto'
        if nnz(interactionMatrix)/numel(interactionMatrix) < 0.5
            interactionMatrix = sparse(interactionMatrix);
        end
    otherwise
        error('Unknown interaction storage: %s.',cfg.numerics.interactionStorage);
end

end
