function [segmentFeatures] = AddEdgeToIntercalate(Geo, numCell, segmentFeatures, edgeLengths_Top, edgesToIntercalate_Top, ghostNodeID)
%ADDEDGETOINTERCALATE Summary of this function goes here
%   Detailed explanation goes here
    for neighbourToNumCell = find(edgesToIntercalate_Top)'
        neighbours_1 = getNodeNeighboursPerDomain(Geo, numCell, ghostNodeID);
        neighbours_2 = getNodeNeighboursPerDomain(Geo, neighbourToNumCell, ghostNodeID);
        sharedNeighbours = intersect(neighbours_1, neighbours_2);

        sharedGhostNodes = sharedNeighbours(ismember(sharedNeighbours, Geo.XgID));

        neighbours_1 = {neighbours_1};
        edgeLength = edgeLengths_Top(neighbourToNumCell);

        for nodePair_g = sharedGhostNodes'
            neighbours_2 = {getNodeNeighboursPerDomain(Geo, nodePair_g, nodePair_g)};
            sharedNeighbours = intersect(neighbours_1{1}, neighbours_2{1});
            sharedNeighbours_c = sharedNeighbours(ismember(sharedNeighbours, Geo.XgID) == 0);
            if any([Geo.Cells(sharedNeighbours_c).AliveStatus] == 1)
                cellToIntercalate = sharedNeighbours_c([Geo.Cells(sharedNeighbours_c).AliveStatus] == 1);
            else
                cellToIntercalate = -1;
            end
            cFace = getFacesFromNode(Geo, [numCell nodePair_g]);
            cFace = cFace{1};
            sharedNeighbours = {sharedNeighbours};
            faceGlobalId = cFace.globalIds;
            segmentFeatures(end+1, :) = table(numCell, nodePair_g, cellToIntercalate, edgeLength, sharedNeighbours, faceGlobalId, neighbours_1, neighbours_2);
        end
    end
end

