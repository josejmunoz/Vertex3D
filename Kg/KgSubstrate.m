function [g,K,Cell,energy] = KgSubstrate(Cell, Y, Yn, Set)
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

        %% Calculate residual g
        g_current = computeGSubstrate(kSubstrate, Y.DataRow(numVertex, 3), Yn.DataRow(numVertex, 3));
        g = Assembleg(g, g_current, numVertex);
        
        %% Save contractile forces (g) to output
        substrateForcesOfCell(numVertexElem, 1) = g_current(3);
        
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKSubstrate(kSubstrate);

            if Set.Sparse
                [si,sj,sv,sk] = AssembleKSparse(K_current, numVertex, si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, numVertex);
            end

            %% Calculate energy
            energy = energy + computeEnergySubstrate(kSubstrate, Y.DataRow(numVertex, 3), Yn.DataRow(numVertex, 3));
        end
    end

    %Cell.SubstrateForce(numCell) = {substrateForcesOfCell};
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

