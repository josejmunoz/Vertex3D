function [g,K,Cell,energy] = KgSubstrate(Cell, SCn, Y, Yn, Set)
%KGSUBSTRATE Summary of this function goes here
%   Detailed explanation goes here
%% Set parameters
Set.Sparse = true;

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
    
    if Cell.GhostCells(numCell)
        kSubstrate = 0;
    else
        kSubstrate = Set.kSubstrate;
    end

    substrateForcesOfCell = zeros(length(Cell.BasalVertices), 1);
    numVertexElem = 1;
    for numVertex = Cell.BasalVertices{numCell}'
        
        if numVertex < 0 %% Face centre
            currentVertex = Cell.FaceCentres.DataRow(abs(numVertex), :);
            currentVertexPrev = SCn.DataRow(abs(numVertex), :);
            vertexIndex = abs(numVertex) + Y.n;
        else %% Regular Vertex
            currentVertex = Y.DataRow(numVertex, :);
            currentVertexPrev = Yn.DataRow(numVertex, :);
            vertexIndex = numVertex;
        end

        %% Calculate residual g
        g_current = computeGSubstrate(kSubstrate, currentVertex(:, 3), currentVertexPrev(:, 3));
        g = Assembleg(g, g_current, vertexIndex);
        
        %% Save contractile forces (g) to output
        substrateForcesOfCell(numVertexElem, 1) = g_current(3);
        numVertexElem = numVertexElem + 1;
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
            energy = energy + computeEnergySubstrate(kSubstrate, currentVertex(:, 3), currentVertexPrev(:, 3));
        end
        
    end

    Cell.SubstrateForce(numCell) = {substrateForcesOfCell};
end

if Set.Sparse && nargout>1
    K = sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
end

end

function [kContractility] = computeKSubstrate(K)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

    kContractility(1:3, 1:3) = [0 0 0; 0 0 0; 0 0 K];

end

function [gContractility] = computeGSubstrate(K, Yz, Yz0)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

gContractility(1:3, 1) = [0 0 (K * (Yz - Yz0))];

end

function [energyConctratility] = computeEnergySubstrate(K, Yz, Yz0)

energyConctratility = 1/2 * K * (Yz - Yz0)^2;

end

