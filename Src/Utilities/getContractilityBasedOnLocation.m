function [contractilityValue] = getContractilityBasedOnLocation(currentFace, currentTri, Geo, Set)
%GETCONTRACTILITYBASEDONLOCATION Summary of this function goes here
%   Detailed explanation goes here
    
    faceConnections = currentFace.ij;
    faceConnections(ismember(faceConnections, Geo.AssembleNodes) == 0) = [];
    timeAfterAblation = Set.currentT - Set.TInitAblation;
    if timeAfterAblation >= 0
        distanceToTimeVariables = abs(Set.Contractility_TimeVariability - timeAfterAblation)/Set.Contractility_TimeVariability(2);
        [closestTimePointsDistance, indicesOfClosestTimePoints] = sort(distanceToTimeVariables);
        closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages
    end
    
    switch (currentFace.InterfaceType)
        case 'Top' % Top
            if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                contractilityValue = Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                    Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);
                contractilityValue = contractilityValue * Set.cLineTension;
            else
                contractilityValue = Set.cLineTension;
            end
        case 'CellCell' % Lateral
            if any([Geo.Cells(currentTri.SharedByCells).AliveStatus] == 0)
                contractilityValue = Set.Contractility_Variability_LateralCables(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                    Set.Contractility_Variability_LateralCables(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);
                contractilityValue = contractilityValue * Set.cLineTension;
            else
                contractilityValue = Set.cLineTension/100;
            end
        case 'Bottom' % Bottom/Substrate
            contractilityValue = Set.cLineTension/100;
        otherwise
            contractilityValue = Set.cLineTension;
    end
end

