function [g,K,Cell,energy] = KgSubstrate(Cell, SCn, Y, Yn, Set)
%KGSUBSTRATE Summary of this function goes here
%   Detailed explanation goes here

%% Initialize
dimensionsG = Set.NumTotalV*3;
g=zeros(dimensionsG,1); % Local cell residual

if Set.Sparse
    sk=0;
    si=zeros((dimensionsG*3)^2,1); % Each vertex is shared by at least 3 cells 
    sj=si;
    sv=si;
    K=sparse(zeros(dimensionsG)); % Also used in sparse
else
    K=zeros(dimensionsG); % Also used in sparse
end

energy = 0;

for numCell = 1:Cell.n
    
    if Cell.DebrisCells(numCell)
        kSubstrate = 0;
        edgeLocation = zeros(size(Cell.EdgeLocation{numCell}));
    else
        edgeLocation = Cell.EdgeLocation{numCell};
        kSubstrate = Set.kSubstrate;
    end

    substrateForcesOfCell = Cell.SubstrateForce{numCell};
    numVertexElem = 0;
    basalJunctionVertices = Cell.BasalBorderVertices{numCell};
    for numVertex = Cell.BasalVertices{numCell}'
        numVertexElem = numVertexElem + 1;
        if basalJunctionVertices(numVertexElem) == 0
            continue;
        end
        
        z0 = Set.z0Substrate;
        if numVertex < 0 %% Face centre
            continue
%             currentVertex = Cell.FaceCentres.DataRow(abs(numVertex), :);
%             %z0 = SCn.DataRow(abs(numVertex), 3);
%             vertexIndex = abs(numVertex) + Set.NumMainV;
%             kSubstrate = 0;
        else %% Regular Vertex
            currentVertex = Y.DataRow(numVertex, :);
            %z0 = Yn.DataRow(numVertex, 3);
            vertexIndex = numVertex;
        end

        %% Calculate residual g
        g_current = computeGSubstrate(kSubstrate, currentVertex(:, 3), z0);
        g = Assembleg(g, g_current, vertexIndex);
        
        %% Save contractile forces (g) to output
        substrateForcesOfCell(numVertexElem, 1) = g_current(3);
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKSubstrate(kSubstrate);

            if Set.Sparse
                [si,sj,sv,sk] = AssembleKSparse(K_current, vertexIndex, si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, vertexIndex);
            end

            %% Calculate energy
            energy = energy + computeEnergySubstrate(kSubstrate, currentVertex(:, 3), z0);
        end
        
    end

    Cell.SubstrateForce(numCell) = {substrateForcesOfCell};
end

if Set.Sparse && nargout>1
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

