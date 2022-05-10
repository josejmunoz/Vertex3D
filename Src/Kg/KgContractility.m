function [g,K,Cell,Energy] = KgContractility(Cell,Y,Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

    %% Initialize
    [g, K] = initializeKg(Geo, Set); 
	
    %% Loop over Cells 
	% Analytical residual g and Jacobian K
	for numCell=1:Geo.nCells
        if Geo.Remodelling
			if ~ismember(numCell,Geo.AssembleNodes)
        		continue
			end
        end
        if Geo.Cells(numCell).AliveStatus ~= 1
            continue
        end
        
    edgeVertices = Cell.Cv{numCell};
    
    edgeLengths = Cell.EdgeLengths{numCell};
    edgeLengths0_average = Cell.EdgeLengths0_average;
    edgeLengths0_lateralAverage = Cell.EdgeLengths0_lateralAverage;
    
    
    edgeLocation = Cell.EdgeLocation{numCell};
    
    contractileForcesOfCell = zeros(size(edgeVertices, 1), 1);
    
    if Set.Sparse > 0
        ge=sparse(size(g, 1), 1); % Local cell residual
    else
        ge=zeros(size(g, 1), 1);
    end
    
    if any(Cell.DebrisCells)
        neighbourWoundEdges = vertcat(Cell.Cv{Cell.DebrisCells});
        %WoundEdge
        idShareEdges = find(ismember(sort(edgeVertices, 2), sort(neighbourWoundEdges, 2), 'rows'));
        
        %FirstRow
        cells1Row = Cell.Int(ismember(Cell.Int, neighbourWoundEdges));
        neighbours1Row = vertcat(Cell.Cv{cells1Row});
        idShare1Row = find(ismember(sort(edgeVertices, 2), sort(neighbours1Row, 2), 'rows'));
        
        %Second Row
        cells2Row = Cell.Int(ismember(Cell.Int, neighbours1Row));
        neighbours2Row = vertcat(Cell.Cv{cells2Row});
        idShare2Row = find(ismember(sort(edgeVertices, 2), sort(neighbours2Row, 2), 'rows'));
        
        %Third Row
        cells3Row = Cell.Int(ismember(Cell.Int, neighbours2Row));
        neighbours3Row = vertcat(Cell.Cv{cells3Row});
        idShare3Row = find(ismember(sort(edgeVertices, 2), sort(neighbours3Row, 2), 'rows'));
    else
        idShareEdges = [];
        idShare1Row = [];
        idShare2Row = [];
        idShare3Row = [];
    end
    
    for numEdge = 1:length(edgeLengths)
        y_1 = Y.DataRow(edgeVertices(numEdge, 1), :);
        l_i0 = edgeLengths0_average;
        if  edgeVertices(numEdge, 2) > 0 %Vertex
            y_2 = Y.DataRow(edgeVertices(numEdge, 2), :);
        else %Face center
            y_2 = Cell.FaceCentres.DataRow(abs(edgeVertices(numEdge, 2)), :);
            edgeVertices(numEdge, 2) = abs(edgeVertices(numEdge, 2)) + Set.NumMainV;
        end
        
        if edgeLocation(numEdge) == 3 % Apical side
            if ismember(numEdge, idShareEdges) % Wound edge
                C = Set.cPurseString;
            else
                C = Set.cLineTension;
            end
        elseif edgeLocation(numEdge) == 2 % Basal side
            C = Set.cLineTension/100;
        elseif edgeLocation(numEdge) == 1  %lateralCables
            edgePosition = ismember(numEdge, idShareEdges) + ismember(numEdge, idShare1Row)...
                 + ismember(numEdge, idShare2Row) + ...
                 ismember(numEdge, idShare3Row);
            switch (edgePosition)
                case 0 %Away from wound
                    C = Set.cLineTension/100;
                case 4 % Wound edge
                    C = Set.cLateralCables;
                case 3 % First row
                    C = Set.cLateralCables*3/4;
                case 2 % Second row
                    C = Set.cLateralCables*2/4;
                case 1 % Third row
                    C = Set.cLateralCables*1/4;
            end
            l_i0 = edgeLengths0_lateralAverage;
        else
            C = 0;
        end
        
        %% Calculate residual g
        g_current = computeGContractility(l_i0, y_1, y_2, C);
        ge = Assembleg(ge, g_current, edgeVertices(numEdge, :));
 
%         K_current = computeKContractility(l_i0, l_i, y_1, y_2, C, Set);
%         delta=1e-6;
%         dim = 3;
%         for i=1:2
%             for j=1:dim
%                 if i == 1
%                     y_1(j) = y_1(j) + delta;
%                 else
%                     y_2(j) = y_2(j) + delta;
%                 end
%                 
%                 gB = computeGContractility(l_i0, l_i, y_1, y_2, C, Set);
%                 col=(i-1)*dim+j;
%                 KB(:,col)=(gB-g_current)/delta;
%                 
%                 
%                 if i == 1
%                     y_1(j) = y_1(j) - delta;
%                 else
%                     y_2(j) = y_2(j) - delta;
%                 end
%             end
%         end
%         if C > -0
%             norm(KB - K_current)
%         end
%         K_current = KB;
        %% Save contractile forces (g) to output
        contractileForcesOfCell(numEdge, 1) = norm(g_current(1:3));
        
        %% AssembleK
        if  nargout>1
            %% Calculate Jacobian
            K_current = computeKContractility(l_i0, y_1, y_2, C);

            if Set.Sparse == 2
                [si,sj,sv,sk] = AssembleKSparse(K_current, edgeVertices(numEdge, :), si, sj, sv, sk);
            else
                K = AssembleK(K, K_current, edgeVertices(numEdge, :));
            end

        end
    end
    
    g = g + ge;

    %% Calculate energy
    Energy = Energy + computeEnergyContractility(l_i0, norm(y_1 - y_2), C);
    
    Cell.ContractileForces{numCell} = contractileForcesOfCell;
end

if Set.Sparse == 2 && nargout>1
    K=sparse(si(1:sk),sj(1:sk),sv(1:sk),dimensionsG,dimensionsG);
end

end

function [kContractility] = computeKContractility(l_i0, y_1, y_2, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

dim = 3;

l_i = norm(y_1 - y_2);

kContractility(1:3, 1:3) = -(C / l_i0) *  (1 / l_i^3 * (y_1 - y_2)' * (y_1 - y_2)) + ((C / l_i0) * eye(dim))/l_i;
kContractility(1:3, 4:6) = -kContractility(1:3, 1:3);
kContractility(4:6, 1:3) = -kContractility(1:3, 1:3);
kContractility(4:6, 4:6) = kContractility(1:3, 1:3);

end

function [gContractility] = computeGContractility(l_i0, y_1, y_2, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

l_i = norm(y_1 - y_2);

gContractility(1:3, 1) = (C / l_i0) * (y_1 - y_2) / l_i;
gContractility(4:6, 1) = -gContractility(1:3);

end

function [energyConctratility] = computeEnergyContractility(l_i0, l_i, C)

energyConctratility = (C / l_i0) * l_i;

end