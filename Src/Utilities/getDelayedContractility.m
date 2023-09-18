function contractilityValue = getDelayedContractility(currentT, purseStringStrength, currentTri, CUTOFF)
        %% THE VALUE OF THE CONTRACTILITY IS THE ONE THAT WAS 6 minutes AGO
        delayMinutes = 6;

        distanceToTimeVariables = (currentT - delayMinutes) - currentTri.EdgeLength_time(:, 1);
        contractilityValue = 0;
        if any(distanceToTimeVariables >= 0)
            [closestTimePointsDistance, indicesOfClosestTimePoints] = sort(abs(distanceToTimeVariables));
            closestTimePointsDistance = 1 - closestTimePointsDistance; %To get percentages
            if sum(closestTimePointsDistance == 1)
                CORRESPONDING_EDGELENGTH_6MINUTES_AGO = currentTri.EdgeLength_time(indicesOfClosestTimePoints(1), 2);
            else
                closestTimePointsDistance = closestTimePointsDistance / sum(closestTimePointsDistance(1:2)); %% Average between the two closest elements
                CORRESPONDING_EDGELENGTH_6MINUTES_AGO = currentTri.EdgeLength_time(indicesOfClosestTimePoints(1), 2) * closestTimePointsDistance(1) + ...
                    currentTri.EdgeLength_time(indicesOfClosestTimePoints(2), 2) * closestTimePointsDistance(2);
            end

            if CORRESPONDING_EDGELENGTH_6MINUTES_AGO <= 0
                CORRESPONDING_EDGELENGTH_6MINUTES_AGO = 0;
            end

            contractilityValue = ((CORRESPONDING_EDGELENGTH_6MINUTES_AGO / currentTri.EdgeLength_time(1, 2)) ^ 4.5) * purseStringStrength;
        end

        if contractilityValue < 1
            contractilityValue = 1;
        end

        % THERE SHOULD BE A CUTTOFF OF MAX OF CONTRACTILITY
        if contractilityValue > CUTOFF || isinf(contractilityValue)
            contractilityValue = CUTOFF;
        end
end