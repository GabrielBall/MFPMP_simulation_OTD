function pCurrent = backward_time_discretisation(p,dt,scheme,localOperator,costVariation,rho,gradientOperator,kernelConvolution)

Nx = numel(p);
identity = speye(Nx);

switch lower(char(scheme))

    case 'imex'
        nonlocal = kernelConvolution(rho.*(gradientOperator*p));
        pPrev = (identity-dt.*localOperator)\(p+dt.*(costVariation-nonlocal));

    case 'backward_euler'
        nonlocal = kernelConvolution(spdiags(rho,0,Nx,Nx)*gradientOperator);
        pPrev = (identity-dt.*localOperator+dt.*nonlocal)\(p+dt.*costVariation);

    otherwise
        error('Unknown adjoint time discretisation scheme: %s.',scheme);

end

end
