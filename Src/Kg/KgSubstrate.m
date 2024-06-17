function [g, K, Energy_T, Geo] = KgSubstrate(Geo, Set)
%KGSUBSTRATE Summary of this function goes here
%   Detailed explanation goes here

    %% Initialize
    [g, K] = initializeKg(Geo, Set); 
    
    Energy_T = 0;
    
    kSubstrate = Set.kSubstrate;
    
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
                if ~isequal(currentFace.InterfaceType, 3)
                    continue
                end
                for currentVertex = unique([currentFace.Tris.Edge currentFace.globalIds])
                    z0 = Set.SubstrateZ;

                    if currentVertex <= length(Geo.Cells(c).globalIds)
                        currentVertexYs = currentCell.Y(currentVertex, :);
                        currentGlobalID = Geo.Cells(c).globalIds(currentVertex);
                    else
                        currentVertexYs = currentFace.Centre;
                        currentGlobalID = currentVertex;
                    end

                    %% Calculate residual g
                    g_current = computeGSubstrate(kSubstrate, currentVertexYs(:, 3), z0);
                    ge = Assembleg(ge, g_current, currentGlobalID);

                    %% Save contractile forces (g) to output
                    Geo.Cells(c).SubstrateG(currentVertex) = g_current(3);

                    %% Calculate Jacobian
                    K_current = computeKSubstrate(kSubstrate);
                    K = AssembleK(K, K_current, currentGlobalID);

                    %% Calculate energy
                    Energy_c = Energy_c + computeEnergySubstrate(kSubstrate, currentVertexYs(:, 3), z0);
                end
            end
            g = g + ge;
            Energy(c) = Energy_c;
        end
    end
    Energy_T = sum(Energy);
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


