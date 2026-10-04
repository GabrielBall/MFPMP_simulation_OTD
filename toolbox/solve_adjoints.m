function [P,Q] = solve_adjoints(y,rho,cfg)

%% Input Gathering

Nt = cfg.numTimeSteps;
Nx = cfg.numCells;
dt = cfg.dt;
dx = cfg.dx;
diffusion = 0.5.*cfg.sigmaX.^2;
boundaryCondition = validatestring(lower(char(cfg.domain.boundary)),{'periodic','zero_flux'});

%% Input Validation

if ~isequal(size(y),[1,Nt+1])
    error('y has the wrong size.');
end

if ~isequal(size(rho),[Nx,Nt+1])
    error('rho has the wrong size.');
end

%% Spatial Operators

gradientOperator = spdiags([-ones(Nx,1),ones(Nx,1)],[-1,1],Nx,Nx)./(2.*dx);
laplacianOperator = spdiags([ones(Nx,1),-2.*ones(Nx,1),ones(Nx,1)],[-1,0,1],Nx,Nx)./dx.^2;

if strcmp(boundaryCondition,'periodic')
    gradientOperator(1,Nx) = -1./(2.*dx);
    gradientOperator(Nx,1) = 1./(2.*dx);
    laplacianOperator(1,Nx) = 1./dx.^2;
    laplacianOperator(Nx,1) = 1./dx.^2;
else
    gradientOperator(1,1) = -1./(2.*dx);
    gradientOperator(Nx,Nx) = 1./(2.*dx);
    laplacianOperator(1,1) = -1./dx.^2;
    laplacianOperator(Nx,Nx) = -1./dx.^2;
end

%% Follower Interaction Convolution

if strcmp(boundaryCondition,'periodic')
    displacement = (0:Nx-1).'.*dx;
    displacement = periodic_displacement(displacement,cfg.domain);
    followerWeight = cfg.followerKernel(abs(displacement),cfg.followerKernelParameters);
    kernelTransform = fft(-followerWeight(:).*displacement);
    fftLength = Nx;
else
    displacement = (-(Nx-1):(Nx-1)).'.*dx;
    followerWeight = cfg.followerKernel(abs(displacement),cfg.followerKernelParameters);
    kernelValues = -followerWeight(:).*displacement;
    fftLength = 2.^nextpow2(3.*Nx-2);
    kernelTransform = fft(kernelValues,fftLength);
end

kernelConvolution = @(field) kernel_convolution(field,kernelTransform,fftLength,Nx,dx,boundaryCondition);

%% Terminal Conditions

P = zeros(Nx,Nt+1);
Q = zeros(1,Nt+1);

[~,terminalGradRho] = cfg.terminalCost(rho(:,end),y(end),cfg,cfg.terminalCostParameters);
P(:,end) = terminalGradRho(:)./dx;
Q(end) = 0;

%% Follower Adjoint

for k = Nt:-1:1
    
    % Drift Calculation
    rhoCurrent = rho(:,k);
    yCurrent = y(k);

    followerVelocity = kernelConvolution(rhoCurrent);

    leaderDisplacement = yCurrent-cfg.cellCentres(:);

    if strcmp(boundaryCondition,'periodic')
        leaderDisplacement = periodic_displacement(leaderDisplacement,cfg.domain);
    end

    leaderWeight = cfg.leaderKernel(abs(leaderDisplacement),cfg.leaderKernelParameters);
    leaderVelocity = cfg.leaderFieldScale.*leaderWeight(:).*leaderDisplacement;
    velocity = followerVelocity+leaderVelocity;
    
    % Upwind Discretisation
    positiveVelocity = max(velocity,0);
    negativeVelocity = min(velocity,0);

    if strcmp(boundaryCondition,'zero_flux')
        negativeVelocity(1) = 0;
        positiveVelocity(end) = 0;
    end

    mainDiagonal = (negativeVelocity-positiveVelocity)./dx;
    upperDiagonal = positiveVelocity(1:end-1)./dx;
    lowerDiagonal = -negativeVelocity(2:end)./dx;

    rows = [(1:Nx).';(1:Nx-1).';(2:Nx).'];
    columns = [(1:Nx).';(2:Nx).';(1:Nx-1).'];
    values = [mainDiagonal;upperDiagonal;lowerDiagonal];

    if strcmp(boundaryCondition,'periodic')
        rows = [rows;1;Nx];
        columns = [columns;Nx;1];
        values = [values;-negativeVelocity(1)./dx;positiveVelocity(end)./dx];
    end
    
    % Total Local Operator ( Non-Local Terms evaluated in
    % backward_time_discertisation depending on scheme)
    advectionOperator = sparse(rows,columns,values,Nx,Nx);
    localOperator = advectionOperator+diffusion.*laplacianOperator;

    [~,runningGradRho] = cfg.runningCost(cfg.time(k),rhoCurrent,yCurrent,cfg,cfg.runningCostParameters);
    costVariation = runningGradRho(:)./dx;

    P(:,k) = backward_time_discretisation(P(:,k+1),dt,'imex',localOperator,costVariation,rhoCurrent,gradientOperator,kernelConvolution);

end

%% Leader Adjoint

for k = Nt:-1:1

    rhoCurrent = rho(:,k);
    yCurrent = y(k);

    leaderDisplacement = yCurrent-cfg.cellCentres(:);

    if strcmp(boundaryCondition,'periodic')
        leaderDisplacement = periodic_displacement(leaderDisplacement,cfg.domain);
    end
    
    % Leader Influence
    leaderDistance = abs(leaderDisplacement);
    [leaderWeight,dLeaderWeightDr] = cfg.leaderKernel(leaderDistance,cfg.leaderKernelParameters);

    % Leader-Induced Velocity Derivative 
    dLeaderVelocityDy = cfg.leaderFieldScale.*(leaderWeight(:)+dLeaderWeightDr(:).*leaderDistance);

    [~,~,runningGradY] = cfg.runningCost(cfg.time(k),rhoCurrent,yCurrent,cfg,cfg.runningCostParameters);

    % Approximate Leader Interaction Integral \int D_y g(y^n)(x) * grad P^n(x) * rho^n(x) dx
    leaderInteraction = dx.*sum(dLeaderVelocityDy.*(gradientOperator*P(:,k)).*rhoCurrent);
    leaderAdjointRate = runningGradY+leaderInteraction;

    % Backward Euler Update
    Q(k) = Q(k+1)+dt.*leaderAdjointRate;

end

end

function convolution = kernel_convolution(field,kernelTransform,fftLength,Nx,dx,boundaryCondition)

field = full(field);

if isvector(field)
    field = field(:);
end

if strcmp(boundaryCondition,'periodic')
    convolution = dx.*real(ifft(fft(field,[],1).*kernelTransform,[],1));
else
    fullConvolution = real(ifft(fft(field,fftLength,1).*kernelTransform,[],1));
    convolution = dx.*fullConvolution(Nx:2.*Nx-1,:);
end

end
