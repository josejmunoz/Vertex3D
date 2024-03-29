function [g, K, Energy_T, Geo] = KgContractility(Geo, Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

    %% Initialize
    [g, K] = initializeKg(Geo, Set);

    %% Reduce dimensionality of K because we won't use any FaceCentre
    % It will increase the efficienciy of 'AssembleK'
    oldSize = size(K, 1);
    K = K(1:Geo.numY*3, 1:Geo.numY*3);

    %% Loop over Cells 
	% Analytical residual g and Jacobian K
	for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        currentCell = Geo.Cells(c);
        if Geo.Remodelling
			if ~ismember(currentCell.ID, Geo.AssembleNodes)
        		continue
			end
        end

        if Geo.Cells(c).AliveStatus
            ge=sparse(size(g, 1), 1);
            Energy_c = 0;
            
            for numFace = 1:length(currentCell.Faces)
                currentFace = Geo.Cells(c).Faces(numFace);
                l_i0 = Geo.EdgeLengthAvg_0(double(currentFace.InterfaceType)+1);
                
                for numTri = 1:length(currentFace.Tris)
                    currentTri = Geo.Cells(c).Faces(numFace).Tris(numTri);
                    if length(currentTri.SharedByCells) > 1
                        [C, Geo] = getContractilityBasedOnLocation(currentFace, currentTri, Geo, Set);
                        
                        y_1 = currentCell.Y(currentTri.Edge(1), :);
                        y_2 = currentCell.Y(currentTri.Edge(2), :);
    
                        %% Calculate residual g
                        g_current = computeGContractility(l_i0, y_1, y_2, C);
                        ge = Assembleg(ge, g_current, currentCell.globalIds(currentTri.Edge));
    
                        %% Save contractile forces (g) to output
                        Geo.Cells(c).Faces(numFace).Tris(numTri).ContractileG = norm(g_current(1:3));
    
                        %% Calculate Jacobian
                        K_current = computeKContractility(l_i0, y_1, y_2, C);
    
                        K = AssembleK(K, K_current, currentCell.globalIds(currentTri.Edge));
    
                        %% Calculate energy
                        Energy_c = Energy_c + computeEnergyContractility(l_i0, norm(y_1 - y_2), C);                    
                    end
                end
            end
            
            g = g + ge;
            Energy(c) = Energy_c;
        end
    end
    
    %% Add zeros and increase size to make it comparable
    K(end+1:oldSize, end+1:oldSize) = 0;

    Energy_T = sum(Energy);
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