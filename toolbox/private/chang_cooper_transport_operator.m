function [A,dLeftDv,dRightDv,dLeaderVelocityDy] = chang_cooper_transport_operator(rho,y,cfg,interactionMatrix,boundaryCondition)

% Velocities
rho = rho(:);
followerVelocity = interactionMatrix*rho;
displacementToLeader = y-cfg.cellEdges(:);

if strcmp(boundaryCondition,'periodic')
    displacementToLeader = periodic_displacement(displacementToLeader,cfg.domain);
end

distanceToLeader = abs(displacementToLeader);

if nargout >= 4
    [leaderWeight,dLeaderWeightDr] = cfg.leaderKernel(distanceToLeader,cfg.leaderKernelParameters);
    dLeaderVelocityDy = cfg.leaderFieldScale.*(leaderWeight+dLeaderWeightDr.*distanceToLeader);
else
    leaderWeight = cfg.leaderKernel(distanceToLeader,cfg.leaderKernelParameters);
    dLeaderVelocityDy = [];
end

leaderVelocity = cfg.leaderFieldScale.*leaderWeight.*displacementToLeader;
interfaceVelocity = followerVelocity+leaderVelocity;

% Flux operator
if nargout >= 3
    [leftCoefficient,rightCoefficient,dLeftDv,dRightDv] = mean_field_flux_coefficients(interfaceVelocity,0.5.*cfg.sigmaX.^2,cfg.dx);
else
    [leftCoefficient,rightCoefficient] = mean_field_flux_coefficients(interfaceVelocity,0.5.*cfg.sigmaX.^2,cfg.dx);
    dLeftDv = [];
    dRightDv = [];
end

A = -build_flux_operator(leftCoefficient,rightCoefficient,cfg.dx,boundaryCondition);

end
