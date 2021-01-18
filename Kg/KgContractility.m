function [g,K,Cell,energy] = KgContractility(Cell,Y,Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

%% Set parameters
C = Set.cContractility;
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

%% Calculate basic information
[Cell] = Cell.computeEdgeLengths(Y);

for numCell = 1:Cell.n
    edgeVertices = Cell.Cv{numCell};
    
    if ~Cell.GhostCells(numCell)
        Cell.ContractileForces{numCell} = zeros(size(edgeVertices, 1), 1);
        continue;
    end
    
    edgeLengths = Cell.EdgeLengths{numCell};
    edgeLengths0_average = Cell.EdgeLengths0_average;
    
    contractileForcesOfCell = zeros(size(edgeVertices, 1), 1);
    
    for numEdge = 1:length(edgeLengths)
        y_1 = Y.DataRow(edgeVertices(numEdge, 1), :);
        
        if  edgeVertices(numEdge, 2) > 0 %Vertex
            y_2 = Y.DataRow(edgeVertices(numEdge, 2), :);
        else %Face center
            y_2 = Cell.FaceCentres.DataRow(abs(edgeVertices(numEdge, 2)), :);
            edgeVertices(numEdge, 2) = abs(edgeVertices(numEdge, 2)) + Y.n;
        end
        
        l_i = edgeLengths(numEdge, 1);
        l_i0 = edgeLengths0_average;
        
        %% Calculate residual g
        g_current = computeGContractility(l_i0, l_i, y_1, y_2, C, Set);
        g = Assembleg(g, g_current, edgeVertices(numEdge, :));
        contractileForcesOfCell(numEdge, 1) = norm(g_current(1:3));
        
        %% Save contractile forces (g) to output
        
        
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKContractility(l_i0, l_i, y_1, y_2, C, Set);

            if Set.Sparse
                [si,sj,sv,sk] = AssembleKSparse(K_current, edgeVertices(numEdge, :), si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, edgeVertices(numEdge, :));
            end

            %% Calculate energy
            energy = energy + computeEnergyContractility(l_i0, l_i, C, Set);
        end
    end
    
    Cell.ContractileForces{numCell} = contractileForcesOfCell;
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
end

end

function [kContractility] = computeKContractility(l_i0, l_i, y_1, y_2, C, Set)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

dim = 3;

if Set.Contractility == 0
    kContractility(1:3, 1:3) = -(C / l_i0) *  ((1 / l_i^2) * (y_1 - y_2)) .* ((1/l_i) * (y_1 - y_2)') + ((C / l_i0) * eye(dim));
    kContractility(1:3, 4:6) = -kContractility(1:3, 1:3);
    kContractility(4:6, 1:3) = -kContractility(1:3, 1:3);
    kContractility(4:6, 4:6) = kContractility(1:3, 1:3);
elseif Set.Contractility == 1
    % Angle added
elseif Set.Contractility == 2
    % Adjacent triangle areas added
end

end

function [gContractility] = computeGContractility(l_i0, l_i, y_1, y_2, C, Set)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

if Set.Contractility == 0
    gContractility(1:3, 1) = (C / l_i0) * (y_1 - y_2) / l_i;
    gContractility(4:6, 1) = -gContractility(1:3);
elseif Set.Contractility == 1
    % Angle added
elseif Set.Contractility == 2
    % Adjacent triangle areas added
end

end

function [energyConctratility] = computeEnergyContractility(l_i0, l_i, C, Set)

if Set.Contractility == 0
    energyConctratility = (C / l_i0) * l_i;
elseif Set.Contractility == 1
    % Angle added
elseif Set.Contractility == 2
    % Adjacent triangle areas added
end

end