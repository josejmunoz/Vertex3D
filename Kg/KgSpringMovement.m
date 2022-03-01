function [g,K,Cell,Energy] = KgSubstrate(Cell, Y, Set)
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
Set.z0Substrate = 4;
for numCell = 1:ncell
    if numCell ~= 10 % Only ductal cell
        continue
    end
    
%     if Cell.DebrisCells(numCell)
%         kSubstrate = 0;
%     else
%         kSubstrate = Set.kSubstrate;
%     end

    if Set.Sparse > 0
        ge=sparse(size(g, 1), 1); % Local cell residual
    else
        ge=zeros(size(g, 1), 1);
    end
    
    kSubstrate = Set.kSubstrate;

    substrateForcesOfCell = Cell.SubstrateForce{numCell};
    numVertexElem = 0;
    %basalJunctionVertices = Cell.BasalBorderVertices{numCell};
    
    currentEdgesOfCell = Cell.Cv{numCell};
    uniqueCurrentVertices = unique(currentEdgesOfCell(currentEdgesOfCell > 0));
    uniqueCurrentFaceCentres = unique(currentEdgesOfCell(currentEdgesOfCell <= 0));
    distances = pdist2(vertcat(Y.DataRow(uniqueCurrentVertices, :), Cell.FaceCentres.DataRow(abs(uniqueCurrentFaceCentres), :)), [0 0 Set.z0Substrate], 'euclidean');
    distances = distances.^10;
    normalizedDistances = (distances)/max(distances);
    normalizedDistances (normalizedDistances ~= 1) = normalizedDistances (normalizedDistances ~= 1).^10;
    for numVertex = vertcat(uniqueCurrentVertices, uniqueCurrentFaceCentres)'
        numVertexElem = numVertexElem + 1;
        currentWeight = normalizedDistances(numVertexElem);
%         if basalJunctionVertices(numVertexElem) == 0
%             continue;
%         end
        
        z0 = Set.z0Substrate;
        if numVertex < 0 %% Face centre
            currentVertex = Cell.FaceCentres.DataRow(abs(numVertex), :);
            vertexIndex = abs(numVertex) + Set.NumMainV;
        else %% Regular Vertex
            currentVertex = Y.DataRow(numVertex, :);
            vertexIndex = numVertex;
        end

        %% Calculate residual g
        g_current = computeGSubstrate(kSubstrate * currentWeight, currentVertex(:, 3), z0);
        ge = Assembleg(ge, g_current, vertexIndex);
        
        %% Save contractile forces (g) to output
        substrateForcesOfCell(numVertexElem, 1) = g_current(3);
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKSubstrate(kSubstrate);

            if Set.Sparse == 2
                [si,sj,sv,sk] = AssembleKSparse(K_current * currentWeight, vertexIndex, si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, vertexIndex);
            end

            %% Calculate energy
            Energy = Energy + computeEnergySubstrate(kSubstrate * currentWeight, currentVertex(:, 3), z0);
        end
    end
    
    g = g + ge;
    Cell.SubstrateForce(numCell) = {substrateForcesOfCell};
end

if Set.Sparse == 2 && nargout>1
    K = sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
end

end

function [kSubstrate] = computeKSubstrate(K)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

    kSubstrate(1:3, 1:3) = [0 0 0; 0 0 0; 0 0 K];

end

function [gSubstrate] = computeGSubstrate(K, Yz, Yz0)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

gSubstrate(1:3, 1) = [0 0 (K * (Yz - Yz0))];

end

function [energySubstrate] = computeEnergySubstrate(K, Yz, Yz0)

energySubstrate = 1/2 * K * (Yz - Yz0)^2;

end

