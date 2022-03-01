function [g,K,Cell,Energy] = KgSpringMovement(Cell, Y, Set)
%KGSUBSTRATE Summary of this function goes here
%   Detailed explanation goes here

%% Initialize
if nargout > 1
    if Set.Sparse == 2 %Manual sparse
        [g, Energy, ncell, K, si, sj, sk, sv] = initializeKg(Cell, Set);
    else %Matlab sparse
        [g, Energy, ncell, K] = initializeKg(Cell, Set);
    end
else
    [g, Energy, ncell] = initializeKg(Cell, Set);
end

z0 = Set.DestinationPoint(3);
for numCell = find(Cell.CellTypes == 2)
    if Set.Sparse > 0
        ge=sparse(size(g, 1), 1); % Local cell residual
    else
        ge=zeros(size(g, 1), 1);
    end
    
    movementStrength = Set.MovementStrength;

    numVertexElem = 0;
    %basalJunctionVertices = Cell.BasalBorderVertices{numCell};
    
    currentEdgesOfCell = Cell.Cv{numCell};
    uniqueCurrentVertices = unique(currentEdgesOfCell(currentEdgesOfCell > 0));
    uniqueCurrentFaceCentres = unique(currentEdgesOfCell(currentEdgesOfCell <= 0));
    distances = pdist2(vertcat(Y.DataRow(uniqueCurrentVertices, :), Cell.FaceCentres.DataRow(abs(uniqueCurrentFaceCentres), :)), Set.DestinationPoint, 'euclidean');
    distances = distances.^10;
    normalizedDistances = (distances)/max(distances);
    normalizedDistances (normalizedDistances ~= 1) = normalizedDistances (normalizedDistances ~= 1).^10;
    for numVertex = vertcat(uniqueCurrentVertices, uniqueCurrentFaceCentres)'
        numVertexElem = numVertexElem + 1;
        currentWeight = normalizedDistances(numVertexElem);
        
        if numVertex < 0 %% Face centre
            currentVertex = Cell.FaceCentres.DataRow(abs(numVertex), :);
            vertexIndex = abs(numVertex) + Set.NumMainV;
        else %% Regular Vertex
            currentVertex = Y.DataRow(numVertex, :);
            vertexIndex = numVertex;
        end

        %% Calculate residual g
        g_current = computeGMovement(movementStrength * currentWeight, currentVertex(:, 3), z0);
        ge = Assembleg(ge, g_current, vertexIndex);
        
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKMovement(movementStrength);

            if Set.Sparse == 2
                [si,sj,sv,sk] = AssembleKSparse(K_current * currentWeight, vertexIndex, si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, vertexIndex);
            end

            %% Calculate energy
            Energy = Energy + computeEnergyMovement(movementStrength * currentWeight, currentVertex(:, 3), z0);
        end
    end
    
    g = g + ge;
end

if Set.Sparse == 2 && nargout>1
    K = sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
end

end

function [kMovement] = computeKMovement(K)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

    kMovement(1:3, 1:3) = [0 0 0; 0 0 0; 0 0 K];

end

function [gMovement] = computeGMovement(K, Yz, Yz0)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

gMovement(1:3, 1) = [0 0 (K * (Yz - Yz0))];

end

function [energyMovement] = computeEnergyMovement(K, Yz, Yz0)

energyMovement = 1/2 * K * (Yz - Yz0)^2;

end

