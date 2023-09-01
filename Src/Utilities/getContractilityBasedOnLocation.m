function [contractilityValue, Geo] = getContractilityBasedOnLocation(currentFace, currentTri, Geo, Set)
%GETCONTRACTILITYBASEDONLOCATION Summary of this function goes here
%   Detailed explanation goes here
    
    noiseContractility = 0.1;
    CUTOFF = 10;
    
    if isempty(currentTri.ContractilityValue)

        %% THE VALUE OF THE CONTRACTILITY IS THE ONE THAT WAS 6 minutes AGO
        delayMinutes = 6;

        distanceToTimeVariables = (Set.currentT - delayMinutes) - currentTri.EdgeLength_time(:, 1);
        contractilityValue = 0;
        if any(distanceToTimeVariables >= 0)
            [closestTimePointsDistance, indicesOfClosestTimePoints] = sort(abs(distanceToTimeVariables));
            closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages
            closestTimePointsDistance = closestTimePointsDistance / sum(closestTimePointsDistance(1:2)); %% Average between the two closest elements
            CORRESPONDING_EDGELENGTH_6MINUTES_AGO = currentTri.EdgeLength_time(indicesOfClosestTimePoints(1), 2) * closestTimePointsDistance(1) + ...
                currentTri.EdgeLength_time(indicesOfClosestTimePoints(2), 2) * closestTimePointsDistance(2);
            contractilityValue = ((CORRESPONDING_EDGELENGTH_6MINUTES_AGO / currentTri.EdgeLength_time(1, 2)) ^ 4.5) * Set.purseStringStrength;
        end

        if contractilityValue < 1
            contractilityValue = 1;
        end

        % THERE SHOULD BE A CUTTOFF OF MAX OF CONTRACTILITY
        if contractilityValue > CUTOFF
            %contractilityValue = CUTOFF;
        end

        switch (currentFace.InterfaceType)
            case 'Top' % Top
                if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                    contractilityValue = contractilityValue * Set.cLineTension;
                else
                    contractilityValue = Set.cLineTension;
                end
            case 'CellCell' % Lateral
                %% DO LATERAL CABLES HAVE A DIFFERENT MINUTES DELAY
                %% CAN IT BE BASED ON HOW FAST IT IS STRAINED?
                if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                    contractilityValue = contractilityValue * Set.cLineTension;
                else
                    contractilityValue = Set.cLineTension/100;
                end
            case 'Bottom' % Bottom/Substrate
                contractilityValue = Set.cLineTension/100;
            otherwise
                contractilityValue = Set.cLineTension;
        end
        %% Adding noise to contractility
        minContractility = contractilityValue - contractilityValue*noiseContractility;
        maxContractility = contractilityValue + contractilityValue*noiseContractility;
        contractilityValue = minContractility + (maxContractility-minContractility)*rand();
    
        for cellToCheck = currentTri.SharedByCells
            facesToCheck = Geo.Cells(cellToCheck).Faces;
            faceToCheckID = ismember(sort(vertcat(facesToCheck.ij), 2), sort(currentFace.ij, 2), "rows");
            if any(faceToCheckID)
                trisToCheck = Geo.Cells(cellToCheck).Faces(faceToCheckID).Tris;
                for triToCheckID = 1:length(trisToCheck)
                    triToCheck = trisToCheck(triToCheckID);
                    if all(ismember(sort(triToCheck.SharedByCells), sort(currentTri.SharedByCells)))
                        Geo.Cells(cellToCheck).Faces(faceToCheckID).Tris(triToCheckID).ContractilityValue = contractilityValue;
                    end
                end
            end
        end
    else
        contractilityValue = currentTri.ContractilityValue;
    end

end

