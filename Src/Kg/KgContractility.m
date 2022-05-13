function [g, K, energy] = KgContractility(Geo, Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

    %% Initialize
    [g, K] = initializeKg(Geo, Set); 
	
    energy = 0;
    %% Loop over Cells 
	% Analytical residual g and Jacobian K
	for numCell = 1:Geo.nCells
        if Geo.Remodelling
			if ~ismember(numCell,Geo.AssembleNodes)
        		continue
			end
        end
        if Geo.Cells(numCell).AliveStatus ~= 1
            continue
        end
        
        ge=sparse(size(g, 1), 1);
        
        currentCell = Geo.Cells(numCell);
        
        [uniqueEdges, uniqueOrder] = unique(vertcat(currentCell.Faces.Tris), 'rows');
        
        edgesLength = vertcat(currentCell.Faces.EdgeLengths);
        edgesLengthOrdered = edgesLength(uniqueOrder);
        
        trisSharedByCells = vertcat(currentCell.Faces.Tris_SharedByCells);
        trisSharedByCells_ordered = trisSharedByCells(uniqueOrder);
        %TODO: GET INFORMATION OF LOCATION OF EDGE (APICAL/BASAL/LATERAL)
        numEdge = 1;
        for currentEdge = uniqueEdges'
            if trisSharedByCells_ordered(numEdge)
                l_i0 = Geo.EdgeLengthsAvg_0;
                y_1 = currentCell.Y(currentEdge(1), :);
                y_2 = currentCell.Y(currentEdge(2), :);

                C = Set.cLineTension;

                %% Calculate residual g
                g_current = computeGContractility(l_i0, y_1, y_2, C);
                ge = Assembleg(ge, g_current, currentCell.globalIds(currentEdge));

                %% Save contractile forces (g) to output
                % TODO
                %contractileForcesOfCell(numEdge, 1) = norm(g_current(1:3));

                %% Calculate Jacobian
                K_current = computeKContractility(l_i0, y_1, y_2, C);

                K = AssembleK(K, K_current, currentCell.globalIds(currentEdge));


                %% Calculate energy
                energy = energy + computeEnergyContractility(l_i0, norm(y_1 - y_2), C);
            end
            numEdge = numEdge + 1;
        end
        
        g = g + ge;

        % TODO:
        %Cell.ContractileForces{numCell} = contractileForcesOfCell;
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