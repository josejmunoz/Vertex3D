function [contractilityValue] = getContractilityBasedOnLocation(currentFace, Geo, Set)
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
        case 0 % Top
            if any([Geo.Cells(faceConnections).AliveStatus] == 0)
                contractilityValue = Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                    Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);
            else
                contractilityValue = Set.cLineTension;
            end
        case 1 % Lateral
            if any([Geo.Cells(faceConnections).AliveStatus] == 0)
                contractilityValue = Set.Contractility_Variability_LateralCables(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                    Set.Contractility_Variability_LateralCables(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);
            else
                contractilityValue = Set.cLineTension;
            end
        case 2 % Bottom/Substrate
            contractilityValue = Set.cLineTension/100;
        otherwise
            contractilityValue = Set.cLineTension;
    end
end
