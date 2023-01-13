function [segmentFeatures] = GetTrisToRemodelOrdered(Geo, Set)
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
for numCell = 1:Geo.nCells
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
                        if cFace.InterfaceType == 'Top'
                            edgeLengths_Top(numSharedCell) = edgeLengths_Top(numSharedCell) + currentTri.EdgeLength;
                        elseif cFace.InterfaceType == 'Bottom'
                            edgeLengths_Bottom(numSharedCell) = edgeLengths_Bottom(numSharedCell) + currentTri.EdgeLength;
                        end
                    end
                end
            end
        end
        
        avgEdgeLengthDomain = mean([Geo.AvgEdgeLength_Bottom, Geo.AvgEdgeLength_Top])*2; %% *2 because there are usually two tris per edge
        edgesToIntercalate_Top = edgeLengths_Top < avgEdgeLengthDomain - (Set.RemodelStiffness * avgEdgeLengthDomain) & edgeLengths_Top > 0;
        edgesToIntercalate_Bottom = edgeLengths_Bottom < avgEdgeLengthDomain - (Set.RemodelStiffness * avgEdgeLengthDomain) & edgeLengths_Bottom > 0;
        
        [segmentFeatures{end+1}] = AddEdgeToIntercalate(Geo, numCell, table(), edgeLengths_Top, edgesToIntercalate_Top, Geo.XgTop(1));
        [segmentFeatures{end+1}] = AddEdgeToIntercalate(Geo, numCell, table(), edgeLengths_Bottom, edgesToIntercalate_Bottom, Geo.XgBottom(1));
    end
end

segmentFeatures(cellfun(@isempty, segmentFeatures)) = [];

end

