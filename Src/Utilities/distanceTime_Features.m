function nonDebris_Features_XMins = distanceTime_Features(Set, timePoints, features_time, XMins)
distanceToTimeVariables = (Set.TInitAblation + XMins) - timePoints;
[closestTimePointsDistance, indicesOfClosestTimePoints] = sort(abs(distanceToTimeVariables));
closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages
if sum(closestTimePointsDistance == 1)
    nonDebris_Features_XMins = features_time{indicesOfClosestTimePoints(1)};
else
    closestTimePointsDistance = closestTimePointsDistance / sum(closestTimePointsDistance(1:2)); %% Average between the two closest elements
    nonDebris_Features_XMins_array = table2array(features_time{indicesOfClosestTimePoints(1)}) * closestTimePointsDistance(1) + ...
        table2array(features_time{indicesOfClosestTimePoints(1)}) * closestTimePointsDistance(2);
    nonDebris_Features_XMins = array2table(nonDebris_Features_XMins_array, 'VariableNames', features_time{1}.Properties.VariableNames);
end
end