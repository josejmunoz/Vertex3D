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
    if ~Cell.GhostCells(numCell)
        continue;
    end
    
    edgeLengths = Cell.EdgeLengths{numCell};
    edgeLengths0_average = Cell.EdgeLengths0_average;
    edgeVertices = Cell.Cv{numCell};
    
    for numEdge = 1:length(edgeLengths)
        y_1 = Y.DataRow(edgeVertices(numEdge, 1), :);
        
        if  edgeVertices(numEdge, 2) > 0 %Vertex
            y_2 = Y.DataRow(edgeVertices(numEdge, 2), :);
        else %Face center
            edgeVertices(numEdge, 2) = abs(edgeVertices(numEdge, 2));
            y_2 = Cell.FaceCentres.DataRow(edgeVertices(numEdge, 2), :);
        end
        
        l_i = edgeLengths(numEdge, 1);
        l_i0 = edgeLengths0_average;
        
        %% Calculate residual g
        g_current = computeGContractility(l_i0, l_i, y_1, y_2, C);
        g = Assembleg(g, g_current, edgeVertices(numEdge, :));
        
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKContractility(l_i0, l_i, y_1, y_2, C);

            if Set.Sparse
                [si,sj,sv,sk] = AssembleKSparse(K_current, edgeVertices(numEdge, :), si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, edgeVertices(numEdge, :));
            end

            %% Calculate energy
            energy = energy + computeEnergyContractility(l_i0, l_i, C);
        end
    end
end

if Set.Sparse && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
end

end

function [kContractility] = computeKContractility(l_i0, l_i, y_1, y_2, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

dim = 3;

%% Check: Is this a matrix?
kContractility(1:3, 1:3) = (C / l_i0) *  ((1 / l_i) * eye(dim) - (1/l_i^3) * (y_1 - y_2) .* (y_1 - y_2)');
kContractility(1:3, 4:6) = -kContractility(1:3, 1:3);
kContractility(4:6, 1:3) = -kContractility(1:3, 1:3);
kContractility(4:6, 4:6) = kContractility(1:3, 1:3);

end

function [gContractility] = computeGContractility(l_i0, l_i, y_1, y_2, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

gContractility(1:3, 1) = (C / l_i0) * (y_1 - y_2) / l_i;
gContractility(4:6, 1) = -gContractility(1:3);

end

function [energyConctratility] = computeEnergyContractility(l_i0, l_i, C)
    energyConctratility = (C / l_i0) * l_i;
end