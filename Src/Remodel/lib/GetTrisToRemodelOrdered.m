function [segmentFeatures] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here

%possibleGhostPairs = nchoosek(Geo.XgID, 2);
allTets = vertcat(Geo.Cells(1:Geo.nCells).T);
allYs = vertcat(Geo.Cells(1:Geo.nCells).Y);

ghostNodeCells = Geo.Cells(ismember(1:length(Geo.Cells), Geo.XgID) & ~cellfun(@isempty, {Geo.Cells.T}));
ghostNodeCellIDs = [ghostNodeCells.ID];
ghostNodeCellIDs = setdiff(ghostNodeCellIDs, Geo.BorderGhostNodes);
ghostPairs = [];
for id = ghostNodeCellIDs
    neighbours = getNodeNeighbours(Geo, id);
    neighbours = neighbours(ismember(neighbours, Geo.XgID));
    idValues = ones(length(neighbours), 1) * id;
    ghostPairs = [ghostPairs; idValues, neighbours];
end

ghostPairs = unique(sort(ghostPairs, 2), 'rows');

segmentFeatures = table();
for ghostPair = ghostPairs'
    ghostPair = ghostPair';
    % Edge length
    x1 = Geo.Cells(ghostPair(1)).X;
    x2 = Geo.Cells(ghostPair(2)).X;
    edgeLength = norm(x2 - x1);
    
    % Number of cell nodes shared
    neighbours_original_1 = getNodeNeighboursPerDomain(Geo, ghostPair(1), ghostPair(1));
    neighbours_original_2 = getNodeNeighboursPerDomain(Geo, ghostPair(2), ghostPair(2));
    neighbours_pair = intersect(neighbours_original_1, neighbours_original_2);
    sharedCellNodes = neighbours_pair(~ismember(neighbours_pair, Geo.XgID));
    neighbours_1 = neighbours_original_1(~ismember(neighbours_original_1, Geo.XgID));
    neighbours_2 = neighbours_original_2(~ismember(neighbours_original_2, Geo.XgID));
    
    neighbours_lengths = [length(neighbours_1), length(neighbours_2)];
    neighbours = {neighbours_1, neighbours_2};
    neighbours_original = {neighbours_original_1, neighbours_original_2};
    
    if any(neighbours_lengths > 3)
        ghostNodesToT1 = ghostPair(neighbours_lengths > 3);
        if length(ghostNodesToT1) > 1
            ghostNodesToT1 = ghostNodesToT1(1);
        end
        neighboursToT1 = neighbours{neighbours_lengths > 3};
        neighboursOriginalToT1 = neighbours_original{neighbours_lengths > 3};
        
        if all(ismember(neighboursToT1, Geo.BorderCells))
            continue
        end
        
        cellNodeNeighbours = {};
        for neigh = neighboursToT1'
            cellNodeNeighbours{end+1} = getNodeNeighboursPerDomain(Geo, ghostNodesToT1, ghostNodesToT1, neigh);
        end
        
        cellNode_Adjacency = cell(length(cellNodeNeighbours), 1);
        for numNode = 1:length(cellNodeNeighbours)
            %otherNodes = setdiff(1:length(cellNodeNeighbours), numNode);
            %opposedNodes{numNode} = setdiff(cellNodeNeighbours{numNode}, vertcat(cellNodeNeighbours{otherNodes}));
            opposedNodes{numNode} = intersect(cellNodeNeighbours{numNode}, Geo.XgID);
            
            currentNeighbours = getNodeNeighboursPerDomain(Geo, neighboursToT1(numNode), ghostNodesToT1);
            cellNode_Adjacency{numNode} = neighboursToT1(ismember(neighboursToT1, [currentNeighbours; neighboursToT1(numNode)]));
        end
        
%         % if one of the neighbours is debris, don't use it to calculate the
%         % mean only
%         cellNode_Adjacency = cellNode_Adjacency([Geo.Cells(neighboursToT1).AliveStatus] > 0);
%         opposedNodes = opposedNodes([Geo.Cells(neighboursToT1).AliveStatus] > 0);
        
        
        if length(cellNode_Adjacency) < 4
            continue
        end
        % Check how they are connected (cell nodes)
        connectedNodes = cellfun(@length, cellNode_Adjacency) == length(neighboursToT1);
        nonConnectedNodes = cellfun(@length, cellNode_Adjacency) < length(neighboursToT1);
        
        nonConnectedXs = neighboursToT1(nonConnectedNodes);
        
        triplet1 = [neighboursToT1(connectedNodes)' nonConnectedXs(1) ghostNodesToT1];
        triplet2 = [neighboursToT1(connectedNodes)' nonConnectedXs(2) ghostNodesToT1];
        
        triplet1_Ys = mean(allYs(sum(ismember(allTets, triplet1), 2) > 3, :));
        triplet2_Ys = mean(allYs(sum(ismember(allTets, triplet2), 2) > 3, :));
        
        distanceEdgeConnectedNodes = norm(triplet1_Ys - triplet2_Ys);
        
%         % Check average distance between connected and unconnected nodes
%         avgDistance_ConnNodes = mean(pdist2(Geo.Cells(ghostNodesToT1).X, vertcat(Geo.Cells(vertcat(opposedNodes{connectedNodes})).X)));
%         avgDistance_NonConnNodes = mean(pdist2(Geo.Cells(ghostNodesToT1).X, vertcat(Geo.Cells(vertcat(opposedNodes{nonConnectedNodes})).X)));
% 
%         if avgDistance_ConnNodes < (avgDistance_NonConnNodes + (Set.RemodelStiffness * avgDistance_NonConnNodes))
%             continue
%         end

        % Check edge length to know when to intercalate
        avgEdgeLengthDomain = mean([Geo.AvgEdgeLength_Bottom, Geo.AvgEdgeLength_Top]);
        if distanceEdgeConnectedNodes > avgEdgeLengthDomain - (Set.RemodelStiffness * avgEdgeLengthDomain)
            continue
        end
        
        % Change the ghostPair to intercalate
        %neighboursToT1 = neighboursToT1([Geo.Cells(neighboursToT1).AliveStatus] > 0);
        if sum(connectedNodes) > 1
            %% T1 transition
            ghostPair = [repmat(ghostNodesToT1, sum(connectedNodes), 1) neighboursToT1(connectedNodes)];
        else
            %% T1 on edge
            ghostPair = [ghostNodesToT1 neighboursToT1(connectedNodes)];
        end
    else
        if edgeLength > Set.RemodelTol
            continue
        end
    end
    
    for gPair = ghostPair'
        % Edge valence (number of shared tets)
        [valence, sharedTets] = edgeValence(Geo, gPair);

        % Face IDs involved
        faceIDs = [];
        for cellNode = sharedCellNodes'
            for numFace = 1:length(Geo.Cells(cellNode).Faces)
                face = Geo.Cells(cellNode).Faces(numFace);
                if any(ismember(gPair, face.ij))
                    faceIDs(end+1) = face.globalIds;
                end
            end
        end

        % Add it to the table
        segmentFeatures(end+1, :) = table(gPair(1), gPair(2), edgeLength, valence, {sharedCellNodes}, {faceIDs}, {neighbours_1}, {neighbours_2}, {sharedTets});
    end
end

if ~isempty(segmentFeatures)
    [segmentFeatures] = sortrows(segmentFeatures, 4, 'ascend');
    [segmentFeatures] = sortrows(segmentFeatures, 3, 'ascend');
end

end

