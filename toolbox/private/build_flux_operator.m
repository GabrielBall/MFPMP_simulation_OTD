function fluxOperator = build_flux_operator( ...
    leftCoefficient,rightCoefficient,dx,boundaryCondition)

% Coefficients
a = leftCoefficient(:);
b = rightCoefficient(:);
numCells = numel(a)-1;
boundaryCondition = lower(char(boundaryCondition));

if numCells == 1
    if strcmp(boundaryCondition,'open')
        fluxOperator = sparse((a(end)+b(1))./dx);
    elseif strcmp(boundaryCondition,'periodic') || strcmp(boundaryCondition,'zero_flux')
        fluxOperator = sparse(1,1);
    else
        error('Unsupported boundary condition: %s.',boundaryCondition);
    end
    return;
end

interiorA = a(2:numCells);
interiorB = b(2:numCells);

% Boundary assembly
switch boundaryCondition
    case {'zero_flux','open'}
        if strcmp(boundaryCondition,'zero_flux')
            leftBoundary = 0;
            rightBoundary = 0;
        else
            leftBoundary = -b(1);
            rightBoundary = a(end);
        end

        mainDiagonal = zeros(numCells,1);
        mainDiagonal(1) = (interiorA(1)-leftBoundary)./dx;
        if numCells > 2
            mainDiagonal(2:end-1) = (interiorA(2:end)+interiorB(1:end-1))./dx;
        end
        mainDiagonal(end) = (rightBoundary+interiorB(end))./dx;

        rows = [(1:numCells).';(1:numCells-1).';(2:numCells).'];
        columns = [(1:numCells).';(2:numCells).';(1:numCells-1).'];
        values = [mainDiagonal;-interiorB./dx;-interiorA./dx];
        fluxOperator = sparse(rows,columns,values,numCells,numCells);

    case 'periodic'
        boundaryA = 0.5.*(a(1)+a(end));
        boundaryB = 0.5.*(b(1)+b(end));

        mainDiagonal = zeros(numCells,1);
        mainDiagonal(1) = (interiorA(1)+boundaryB)./dx;
        if numCells > 2
            mainDiagonal(2:end-1) = (interiorA(2:end)+interiorB(1:end-1))./dx;
        end
        mainDiagonal(end) = (boundaryA+interiorB(end))./dx;

        rows = [(1:numCells).';(1:numCells-1).';(2:numCells).';1;numCells];
        columns = [(1:numCells).';(2:numCells).';(1:numCells-1).';numCells;1];
        values = [mainDiagonal;-interiorB./dx;-interiorA./dx;-boundaryA./dx;-boundaryB./dx];
        fluxOperator = sparse(rows,columns,values,numCells,numCells);

    otherwise
        error('Unsupported boundary condition: %s.',boundaryCondition);
end

end
