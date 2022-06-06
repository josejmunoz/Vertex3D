function [g, K, energy, Geo] = KgContractility(Geo, Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

    %% Initialize
    [g, K] = initializeKg(Geo, Set); 
	
    energy = 0;
    %% Loop over Cells 
	% Analytical residual g and Jacobian K
	for numCell = [Geo.Cells.ID]
        currentCell = Geo.Cells(numCell);
        if Geo.Remodelling
			if ~ismember(currentCell.ID, Geo.AssembleNodes)
        		continue
			end
        end
        if isempty(currentCell.AliveStatus) || currentCell.AliveStatus ~= 1
            continue
        end
        
        ge=sparse(size(g, 1), 1);
        
        for numFace = 1:length(currentCell.Faces)
            currentFace = Geo.Cells(numCell).Faces(numFace);
            l_i0 = Geo.EdgeLengthAvg_0(double(currentFace.InterfaceType)+1);
            
            for numTri = 1:length(currentFace.Tris)
                currentTri = Geo.Cells(numCell).Faces(numFace).Tris(numTri);
                if length(currentTri.SharedByCells) > 1
                    C = getContractilityBasedOnLocation(currentFace, currentTri, Geo, Set);
                    
                    y_1 = currentCell.Y(currentTri.Edge(1), :);
                    y_2 = currentCell.Y(currentTri.Edge(2), :);

                    %% Calculate residual g
                    g_current = computeGContractility(l_i0, y_1, y_2, C);
                    ge = Assembleg(ge, g_current, currentCell.globalIds(currentTri.Edge));

                    %% Save contractile forces (g) to output
                    Geo.Cells(numCell).Faces(numFace).Tris(numTri).ContractileG = norm(g_current(1:3));

                    %% Calculate Jacobian
                    K_current = computeKContractility(l_i0, y_1, y_2, C);

                    K = AssembleK(K, K_current, currentCell.globalIds(currentTri.Edge));

                    %% Calculate energy
                    energy = energy + computeEnergyContractility(l_i0, norm(y_1 - y_2), C);                    
                end
            end
        end
        
        g = g + ge;
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