function [g, K, energy, Geo] = KgSubstrate(Geo, Set)
%KGSUBSTRATE Summary of this function goes here
%   Detailed explanation goes here

    %% Initialize
    [g, K] = initializeKg(Geo, Set); 
    
    energy = 0;
    
    kSubstrate = Set.kSubstrate;
    
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
            if ~isequal(currentFace.InterfaceType, 'Bottom')
                continue
            end
            for currentVertex = unique([currentFace.Tris.Edge])
                z0 = Set.SubstrateZ;
                currentVertexYs = currentCell.Y(currentVertex, :);

                %% Calculate residual g
                g_current = computeGSubstrate(kSubstrate, currentVertexYs(:, 3), z0);
                ge = Assembleg(ge, g_current, Geo.Cells(numCell).globalIds(currentVertex));

                %% Save contractile forces (g) to output
                Geo.Cells(numCell).SubstrateG(currentVertex) = g_current(3);

                %% Calculate Jacobian
                K_current = computeKSubstrate(kSubstrate);
                K = AssembleK(K, K_current, Geo.Cells(numCell).globalIds(currentVertex));

                %% Calculate energy
                energy = energy + computeEnergySubstrate(kSubstrate, currentVertexYs(:, 3), z0);
            end
        end
        g = g + ge;
    end
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


