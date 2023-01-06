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
for id = ghostNodeCellIDs
    neighbours = getNodeNeighbours(Geo, id);
    neighbours = neighbours(ismember(neighbours, ghostNodesWithoutDebris));
    idValues = ones(length(neighbours), 1) * id;
    ghostPairs = [ghostPairs; idValues, neighbours];
end

ghostPairs = unique(sort(ghostPairs, 2), 'rows');

segmentFeatures = table();
for numCell = 1:Geo.nCells
    if Geo.Cells(numCell).AliveStatus
        neighbours_Top = getNodeNeighboursPerDomain(Geo, numCell, Geo.XgTop(1));
        neighbours_Bottom = getNodeNeighboursPerDomain(Geo, numCell, Geo.XgBottom(1));

        Ys = Geo.Cells(numCell).Y;
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
        
        if any(edgesToIntercalate_Top)
            for numCellToIntercalate = find(edgesToIntercalate_Top)'
                neighbours_1 = getNodeNeighboursPerDomain(Geo, numCell, Geo.XgTop(1));
                neighbours_2 = getNodeNeighboursPerDomain(Geo, numCellToIntercalate, Geo.XgTop(1));
                sharedNeighbours = intersect(neighbours_1, neighbours_2);
                
                sharedCellNodes = sharedNeighbours(~ismember(sharedNeighbours, Geo.XgID));
                
                segmentFeatures(end+1, :) = table(numCell, numCellToIntercalate, edgeLengths_Top(edgesToIntercalate_Top), {sharedNeighbours}, {neighbours_1}, {neighbours_2});
            end
        end
        
        if any(edgesToIntercalate_Bottom)
            for numCellToIntercalate = find(edgesToIntercalate_Bottom)'
                neighbours_1 = getNodeNeighboursPerDomain(Geo, numCell, Geo.XgBottom(1));
                neighbours_2 = getNodeNeighboursPerDomain(Geo, numCellToIntercalate, Geo.XgBottom(1));
                sharedNeighbours = intersect(neighbours_1, neighbours_2);
                
                sharedCellNodes = sharedNeighbours(~ismember(sharedNeighbours, Geo.XgID));
                segmentFeatures(end+1, :) = table(numCell, numCellToIntercalate, edgeLengths_Top(edgesToIntercalate_Top), {sharedNeighbours}, {neighbours_1}, {neighbours_2});
            end
        end
    end
end


% for ghostPair = ghostPairs'
%     ghostPair = ghostPair';
%     if any(ismember(ghostPair, Geo.BorderGhostNodes))
%         continue
%     end
%     % Edge length
%     x1 = Geo.Cells(ghostPair(1)).X;
%     x2 = Geo.Cells(ghostPair(2)).X;
%     edgeLength = norm(x2 - x1);
%     
%     % Number of cell nodes shared
%     neighbours_original_1 = getNodeNeighboursPerDomain(Geo, ghostPair(1), ghostPair(1));
%     neighbours_original_2 = getNodeNeighboursPerDomain(Geo, ghostPair(2), ghostPair(2));
%     
%     if any(ismember(neighbours_original_1, Geo.BorderCells)) || any(ismember(neighbours_original_2, Geo.BorderCells))
%         continue
%     end
%     
%     neighbours_pair = intersect(neighbours_original_1, neighbours_original_2);
%     sharedCellNodes = neighbours_pair(~ismember(neighbours_pair, ghostNodesWithoutDebris));
%     neighbours_1 = neighbours_original_1(~ismember(neighbours_original_1, ghostNodesWithoutDebris));
%     neighbours_2 = neighbours_original_2(~ismember(neighbours_original_2, ghostNodesWithoutDebris));
%     
%     neighbours_lengths = [length(neighbours_1), length(neighbours_2)];
%     neighbours = {neighbours_1, neighbours_2};
%     neighbours_original = {neighbours_original_1, neighbours_original_2};
%     
% %     if any(neighbours_lengths > 3)
% %         ghostNodesToT1 = ghostPair(neighbours_lengths > 3);
% %         if length(ghostNodesToT1) > 1
% %             ghostNodesToT1 = ghostNodesToT1(1);
% %         end
% %         neighboursToT1 = neighbours{neighbours_lengths > 3};
% %         neighboursOriginalToT1 = neighbours_original{neighbours_lengths > 3};
% %         
% %         if all(ismember(neighboursToT1, Geo.BorderCells))
% %             continue
% %         end
% %         
% %         aliveCells = ~cellfun(@isempty, {Geo.Cells(neighboursToT1).AliveStatus});
% %         aliveCells(aliveCells) = [Geo.Cells(neighboursToT1).AliveStatus] > 0;
% %         
% %         if sum(aliveCells) < 3
% %             continue
% %         end
% %         
% %         cellNodeNeighbours = {};
% %         for neigh = neighboursToT1'
% %             cellNodeNeighbours{end+1} = getNodeNeighboursPerDomain(Geo, ghostNodesToT1, ghostNodesToT1, neigh);
% %         end
% %         
% %         cellNode_Adjacency = cell(length(cellNodeNeighbours), 1);
% %         for numNode = 1:length(cellNodeNeighbours)
% %             %otherNodes = setdiff(1:length(cellNodeNeighbours), numNode);
% %             %opposedNodes{numNode} = setdiff(cellNodeNeighbours{numNode}, vertcat(cellNodeNeighbours{otherNodes}));
% %             opposedNodes{numNode} = intersect(cellNodeNeighbours{numNode}, ghostNodesWithoutDebris);
% %             
% %             currentNeighbours = getNodeNeighboursPerDomain(Geo, neighboursToT1(numNode), ghostNodesToT1);
% %             cellNode_Adjacency{numNode} = neighboursToT1(ismember(neighboursToT1, [currentNeighbours; neighboursToT1(numNode)]));
% %         end
% %         
% % %         % if one of the neighbours is debris, don't use it to calculate the
% % %         % mean only
% % %         cellNode_Adjacency = cellNode_Adjacency([Geo.Cells(neighboursToT1).AliveStatus] > 0);
% % %         opposedNodes = opposedNodes([Geo.Cells(neighboursToT1).AliveStatus] > 0);
% %         
% %         
% %         if length(cellNode_Adjacency) < 4
% %             continue
% %         end
% %         % Check how they are connected (cell nodes)
% %         connectedNodes = cellfun(@length, cellNode_Adjacency) == length(neighboursToT1);
% %         nonConnectedNodes = cellfun(@length, cellNode_Adjacency) < length(neighboursToT1);
% %         
% %         nonConnectedXs = neighboursToT1(nonConnectedNodes);
% %         
% %         triplet1 = [neighboursToT1(connectedNodes)' nonConnectedXs(1) ghostNodesToT1];
% %         triplet2 = [neighboursToT1(connectedNodes)' nonConnectedXs(2) ghostNodesToT1];
% %         
% %         triplet1_Ys = mean(allYs(sum(ismember(allTets, triplet1), 2) > 3, :));
% %         triplet2_Ys = mean(allYs(sum(ismember(allTets, triplet2), 2) > 3, :));
% %         
% %         if any(isnan(triplet1_Ys)) || any(isnan(triplet2_Ys))
% %             continue
% %         end
% %         
% %         distanceEdgeConnectedNodes = norm(triplet1_Ys - triplet2_Ys);
% % 
% %         % Check edge length to know when to intercalate
% %         avgEdgeLengthDomain = mean([Geo.AvgEdgeLength_Bottom, Geo.AvgEdgeLength_Top]);
% %         if distanceEdgeConnectedNodes > avgEdgeLengthDomain - (Set.RemodelStiffness * avgEdgeLengthDomain)
% %             continue
% %         end
% %         
% %         % Change the ghostPair to intercalate
% %         if sum(connectedNodes) > 1
% %             %% T1 transition
% %             ghostPair = [repmat(ghostNodesToT1, sum(connectedNodes & aliveCells'), 1) neighboursToT1(connectedNodes & aliveCells')];
% %         else
% %             %% T1 on edge
% %             ghostPair = [ghostNodesToT1 neighboursToT1(connectedNodes & aliveCells')];
% %         end
% %     else
%         minEdgeLengths = [];
%         for numGhost = 1:length(ghostPair)
%             currentFaces = getFacesFromNode(Geo, ghostPair(numGhost));
%             neighboursOtherPair = getNodeNeighbours(Geo, setdiff(ghostPair, ghostPair(numGhost)));
%             for f = 1:length(currentFaces)
%                 currentFace = currentFaces{f};
%                 mainNode = setdiff(currentFace.ij, ghostPair(numGhost));
%                 if ismember(mainNode, neighboursOtherPair)
%                     Ys = Geo.Cells(mainNode).Y;
%                     [edgeLengths, LengthsToCentre, AspectRatio, sharedByCells] = ComputeFaceEdgeLengths(currentFace, Ys, 1);
%                     edgeLengths(cellfun(@length, sharedByCells) == 1) = []; 
%                     if length(edgeLengths) > 0
%                         minEdgeLengths(end+1) = min([edgeLengths{:}]);
%                     end
%                 end
%             end
%         end
%             
%         avgEdgeLengthDomain = mean([Geo.AvgEdgeLength_Bottom, Geo.AvgEdgeLength_Top]);
%         if min(minEdgeLengths) > avgEdgeLengthDomain - (Set.RemodelStiffness * avgEdgeLengthDomain)
%             if edgeLength > Set.RemodelTol
%                 continue
%             end
%         end
% %     end
%    
%     %% TODO: AVOID GHOST NODES INTERCALATIONS
%     
%     for gPair = ghostPair'
%         % Edge valence (number of shared tets)
%         [valence, sharedTets] = edgeValence(Geo, gPair);
% 
%         % Face IDs involved
%         faceIDs = [];
%         for cellNode = sharedCellNodes'
%             for numFace = 1:length(Geo.Cells(cellNode).Faces)
%                 face = Geo.Cells(cellNode).Faces(numFace);
%                 if any(ismember(gPair, face.ij))
%                     faceIDs(end+1) = face.globalIds;
%                 end
%             end
%         end
% 
%         % Add it to the table
%         segmentFeatures(end+1, :) = table(gPair(1), gPair(2), edgeLength, valence, {sharedCellNodes}, {faceIDs}, {neighbours_1}, {neighbours_2}, {sharedTets});
%     end
% end

if ~isempty(segmentFeatures)
    %[segmentFeatures] = sortrows(segmentFeatures, 4, 'ascend');
    [segmentFeatures] = sortrows(segmentFeatures, 3, 'ascend');
end

end

