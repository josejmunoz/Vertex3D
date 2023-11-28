function [segmentFeatures_filtered] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here

%possibleGhostPairs = nchoosek(Geo.XgID, 2);
nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
allTets = vertcat(Geo.Cells(nonDeadCells).T);
allYs = vertcat(Geo.Cells(nonDeadCells).Y);
ghostNodesWithoutDebris = setdiff(Geo.XgID, Geo.RemovedDebrisCells);

ghostNodeCells = Geo.Cells(ismember(1:length(Geo.Cells), ghostNodesWithoutDebris) & ~cellfun(@isempty, {Geo.Cells.T}));
ghostNodeCellIDs = [ghostNodeCells.ID];
ghostNodeCellIDs = setdiff(ghostNodeCellIDs, Geo.BorderGhostNodes);
ghostPairs = [];

segmentFeatures = {};
for numCell = nonDeadCells
    if Geo.Cells(numCell).AliveStatus && ~ismember(numCell, Geo.BorderCells)
        
        %%%
        currentFaces = getFacesFromNode(Geo, numCell);
        edgeLengths_Top = zeros(Geo.nCells, 1);
        edgeLengths_Bottom = zeros(Geo.nCells, 1);
        for currentFace = currentFaces
            cFace = currentFace{1};
            currentTris = cFace.Tris;
            for currentTri = currentTris
                if length(currentTri.SharedByCells) > 1
                    sharedCells = currentTri.SharedByCells;
                    sharedCells(sharedCells == numCell) = [];
                    for numSharedCell = sharedCells
                        if cFace.InterfaceType == 1
                            edgeLengths_Top(numSharedCell) = edgeLengths_Top(numSharedCell) + currentTri.EdgeLength / cFace.Area;
                        elseif cFace.InterfaceType == 3
                            edgeLengths_Bottom(numSharedCell) = edgeLengths_Bottom(numSharedCell) + currentTri.EdgeLength / cFace.Area;
                        end
                    end
                end
            end 
        end

        if any(edgeLengths_Top>0)
            avgEdgeLength = mean(edgeLengths_Top(edgeLengths_Top>0));
            edgesToIntercalate_Top = edgeLengths_Top < avgEdgeLength - (Set.RemodelStiffness * avgEdgeLength) & edgeLengths_Top > 0;
            [segmentFeatures{end+1}] = AddEdgeToIntercalate(Geo, numCell, table(), edgeLengths_Top, edgesToIntercalate_Top, Geo.XgTop(1));
        end

        if any(edgeLengths_Bottom>0)
            %avgEdgeLength = Geo.AvgEdgeLength_Bottom/lastFaceArea_Bottom;
            avgEdgeLength = mean(edgeLengths_Bottom(edgeLengths_Bottom>0));
            edgesToIntercalate_Bottom = edgeLengths_Bottom < avgEdgeLength - (Set.RemodelStiffness * avgEdgeLength) & edgeLengths_Bottom > 0;
            [segmentFeatures{end+1}] = AddEdgeToIntercalate(Geo, numCell, table(), edgeLengths_Bottom, edgesToIntercalate_Bottom, Geo.XgBottom(1));
        end
    end
end

segmentFeatures(cellfun(@isempty, segmentFeatures)) = [];
segmentFeatures_filtered = {};
for segmentFeature = segmentFeatures
    segmentFeature = segmentFeature{1};
    gNodeNeighbours = {};
    for numRow = 1:size(segmentFeature, 1)
        gNodeNeighbours{numRow} = getNodeNeighbours(Geo, segmentFeature{numRow, 2});
    end
    gNodes_NeighboursShared = unique(vertcat(gNodeNeighbours{:}));
    cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
    if sum([Geo.Cells(cellNodesShared).AliveStatus]==0) < 2 && length(cellNodesShared) > 3 && length(unique(segmentFeature.cellToSplitFrom)) == 1
        segmentFeatures_filtered{end+1} = segmentFeature;
    end
end
end

