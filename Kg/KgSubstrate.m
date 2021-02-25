function [g,K,Cell,energy] = KgSubstrate(Cell, Y, Set)
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
kSubstrate = Set.kSubstrate;

for numCell = 1:Cell.n

    substrateForce = zeros(length(Cell.BasalVertices), 1);
    
    for numVertex = Cell.BasalVertices(numCell)

        %% Calculate residual g
        g_current = computeGSubstrate(kSubstrate, Y, Yz0);
        g = Assembleg(g, g_current, numVertex);
        
        %% Save contractile forces (g) to output
        substrateForcesOfCell(numVertex, 1) = g_current(3);
        
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
            energy = energy + computeEnergySubstrate(kSubstrate, Yz, Yz0);
        end
    end

    %Cell.ContractileForces{numCell} = contractileForcesOfCell;
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

