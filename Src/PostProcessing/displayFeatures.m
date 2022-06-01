function [measurementsToDisplay_Header, measurementsToDisplay] = displayFeatures(Geo, features, features0, c, featuresToDisplay)
%DISPLAYFEATURES Summary of this function goes here
%   Detailed explanation goes here
    measurementsToDisplay_Header = struct();
    measurementsToDisplay{c} = struct();
    
    for feature = featuresToDisplay'
        numTris = 1;
        measurementsToDisplay_Header.(feature{1}) = "SCALARS " + feature{1} + " double\n";
        measurementsToDisplay_Header.(feature{1}) = measurementsToDisplay_Header.(feature{1}) + "LOOKUP_TABLE default\n";

        for f = 1:length(Geo.Cells(c).Faces)
            %Divide by location when required
            if endsWith(feature{1}, 'ByLocation')
                baseFeature = strrep(feature{1}, 'ByLocation', '');
                currentFeature = strcat(baseFeature, '_', string(Geo.Cells(c).Faces(f).InterfaceType));

                if all(~contains(featuresToDisplay, currentFeature))
                    currentFeature = baseFeature;
                end
            else % Otherwise, we use the same feature
                currentFeature = feature{1};
            end
            for t = 1:length(Geo.Cells(c).Faces(f).Tris)
                %Print the values of the feature regarding the
                %triangle/edge
                if contains(feature{1}, "Tilting") || isempty(features0) % Contains 0s at the beggining and then we get Inf and NaNs
                    result = features(numTris).(currentFeature);
                else
                    result = (features(numTris).(currentFeature) - features0(numTris).(currentFeature)) / features0(numTris).(currentFeature);
                end

                if isfield(measurementsToDisplay{c}, feature{1})
                    measurementsToDisplay{c}.(feature{1}) = measurementsToDisplay{c}.(feature{1}) + sprintf("%f\n", result);
                else
                    measurementsToDisplay{c}.(feature{1}) = sprintf("%f\n", result);
                end
                
                numTris = numTris + 1;
            end
        end
    end
end

