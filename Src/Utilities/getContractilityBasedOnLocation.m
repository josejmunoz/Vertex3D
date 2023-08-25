function [contractilityValue, Geo] = getContractilityBasedOnLocation(currentFace, currentTri, Geo, Set)
%GETCONTRACTILITYBASEDONLOCATION Summary of this function goes here
%   Detailed explanation goes here
    
    noiseContractility = 0.1;
    
    if isempty(currentTri.ContractilityValue)
        timeAfterAblation = Set.currentT - Set.TInitAblation;
        if timeAfterAblation >= 0
            distanceToTimeVariables = abs(Set.Contractility_TimeVariability - timeAfterAblation)/Set.Contractility_TimeVariability(2);
            [closestTimePointsDistance, indicesOfClosestTimePoints] = sort(distanceToTimeVariables);
            closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages
            contractilityValue = Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                        Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);

            %% THE VALUE OF THE CONTRACTILITY IS THE ONE THAT WAS 6 SECONDS AGO
            %% USE THE SAME FUNCTION IN CONTRACTILITY TO OBTAIN THE 
            CORRESPONDING_EDGELENGTH_6SECONDS_AGO = 
            purseString_theory = (CORRESPONDING_EDGELENGTH_6SECONDS_AGO / currentTri.EdgeLength_time(1, 2)) ^ 4.5;
            % THERE SHOULD BE A CUTTOFF OF MAX OF CONTRACTILITY
        end

        switch (currentFace.InterfaceType)
            case 'Top' % Top
                if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                    


                    contractilityValue = contractilityValue * Set.cLineTension;
                else
                    contractilityValue = Set.cLineTension;
                end
            case 'CellCell' % Lateral
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

