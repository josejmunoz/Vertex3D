function [heightLateral, distanceFacesTopBottom] = ComputeCellHeight(Cell)
%COMPUTECELLHEIGHT Summary of this function goes here
%   Detailed explanation goes here
    
    %% Get the height as an average of the distance between the faces of a
    % cell
    allTopFaceCentres = [];
    allBottomFaceCentres = [];
    allLateralFaceCentres = [];
    for face = Cell.Faces
        if face.InterfaceType == 'Top'
            allTopFaceCentres = [allTopFaceCentres; face.Centre];
        elseif face.InterfaceType == 'Bottom'
            allBottomFaceCentres = [allBottomFaceCentres; face.Centre];
        else
            allLateralFaceCentres = [allLateralFaceCentres; face];
        end
    end
    distanceFacesTopBottom = pdist2(mean(allTopFaceCentres), mean(allBottomFaceCentres));
    
    %% Get the height as the length of the lateral edges
    lateralEdgesLength = [];
    for lateralFace = allLateralFaceCentres'
        for tris = lateralFace.Tris
            if length(tris.SharedByCells) > 2
                lateralEdgesLength = [lateralEdgesLength; tris.EdgeLength];
            end
        end
    end
    heightLateral = mean(lateralEdgesLength);
end
