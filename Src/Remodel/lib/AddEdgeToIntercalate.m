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
            sharedNeighbours_c(sharedNeighbours_c == neighbourToNumCell) = [];
            if any([Geo.Cells(sharedNeighbours_c).AliveStatus] == 1)
                cellToIntercalate = sharedNeighbours_c([Geo.Cells(sharedNeighbours_c).AliveStatus] == 1);
            else
                cellToIntercalate = -1;
                continue
            end
            cFace = getFacesFromNode(Geo, [numCell nodePair_g]);
            cFace = cFace{1};
            numSharedNeighbours = length(sharedNeighbours);
            sharedNeighbours = {sharedNeighbours};
            faceGlobalId = cFace.globalIds;
            cellToSplitFrom = neighbourToNumCell;
            for cell_intercalate = cellToIntercalate'
                segmentFeatures(end+1, :) = table(numCell, nodePair_g, cell_intercalate, cellToSplitFrom, edgeLength, numSharedNeighbours, sharedNeighbours, faceGlobalId, neighbours_1, neighbours_2);
            end
%             %% Now numCell should receive a tet from cellToIntercalate
%             neighbours_2 = {getNodeNeighboursPerDomain(Geo, cellToIntercalate, nodePair_g)};
%             sharedNeighbours = intersect(sharedNeighbours{1}, neighbours_2{1});
%             nodeToTransfer = sharedNeighbours(ismember(sharedNeighbours, Geo.XgID));
%             if length(nodeToTransfer)>1
%                 error('culo');
%             end
%             sharedNeighbours = {sharedNeighbours};
%             segmentFeatures(end+1, :) = table(cellToIntercalate, nodeToTransfer, numCell, -1, edgeLength, sharedNeighbours, faceGlobalId, neighbours_1, neighbours_2);
        end
    end
end

