function contractilityValue = getIntensityBasedContractility(Set, currentFace)
    timeAfterAblation = Set.currentT - Set.TInitAblation;
    if timeAfterAblation >= 0
        distanceToTimeVariables = abs(Set.Contractility_TimeVariability - timeAfterAblation)/Set.Contractility_TimeVariability(2);
        [closestTimePointsDistance, indicesOfClosestTimePoints] = sort(distanceToTimeVariables);
        closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages

        switch (currentFace.InterfaceType)
            case 'Top' % Top
                contractilityValue = Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                    Set.Contractility_Variability_PurseString(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);
            case 'CellCell'
                contractilityValue = Set.Contractility_Variability_LateralCables(indicesOfClosestTimePoints(1)) * closestTimePointsDistance(1) + ...
                    Set.Contractility_Variability_LateralCables(indicesOfClosestTimePoints(2)) * closestTimePointsDistance(2);
        end
    end
end