function [segmentFeatures] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here

%possibleGhostPairs = nchoosek(Geo.XgID, 2);

ghostNodeCells = Geo.Cells(ismember(1:length(Geo.Cells), Geo.XgID) & ~cellfun(@isempty, {Geo.Cells.T}));
ghostNodeCellIDs = [ghostNodeCells.ID];
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
    % Edge length
    x1 = Geo.Cells(ghostPair(1)).X;
    x2 = Geo.Cells(ghostPair(2)).X;
    edgeLength = norm(x2 - x1);
    
    % Number of cell nodes shared
    neighbours_original_1 = getNodeNeighboursPerDomain(Geo, ghostPair(1), ismember(ghostPair(1), Geo.XgBottom));
    neighbours_original_2 = getNodeNeighboursPerDomain(Geo, ghostPair(2), ismember(ghostPair(2), Geo.XgBottom));
    neighbours_pair = intersect(neighbours_original_1, neighbours_original_2);
    sharedCellNodes = neighbours_pair(~ismember(neighbours_pair, Geo.XgID));
    neighbours_1 = neighbours_original_1(~ismember(neighbours_original_1, Geo.XgID));
    neighbours_2 = neighbours_original_2(~ismember(neighbours_original_2, Geo.XgID));
    
    neighbours_lengths = [length(neighbours_1), length(neighbours_2)];
    neighbours = {neighbours_1, neighbours_2};
    neighbours_original = {neighbours_original_1, neighbours_original_2};
    
    if any(neighbours_lengths > 3)
        ghostNodesToT1 = ghostPair(neighbours_lengths > 3);
        neighboursToT1 = neighbours{neighbours_lengths > 3};
        neighboursOriginalToT1 = neighbours_original{neighbours_lengths > 3};
        
        if all(ismember(neighboursToT1, Geo.BorderCells))
            continue
        end
        
        cellNodeNeighbours = {};
        for neigh = neighboursToT1'
            cellNodeNeighbours{end+1} = getNodeNeighboursPerDomain(Geo, ghostNodesToT1, ismember(ghostNodesToT1, Geo.XgBottom), neigh);
        end
        
        cellNode_Adjacency = cell(length(cellNodeNeighbours), 1);
        for numNode = 1:length(cellNodeNeighbours)
            %otherNodes = setdiff(1:length(cellNodeNeighbours), numNode);
            %opposedNodes{numNode} = setdiff(cellNodeNeighbours{numNode}, vertcat(cellNodeNeighbours{otherNodes}));
            opposedNodes{numNode} = intersect(cellNodeNeighbours{numNode}, Geo.XgID);
            
            currentNeighbours = getNodeNeighboursPerDomain(Geo, neighboursToT1(numNode), ismember(ghostNodesToT1, Geo.XgBottom));
            cellNode_Adjacency{numNode} = neighboursToT1(ismember(neighboursToT1, [currentNeighbours; neighboursToT1(numNode)]));
        end
        
        % if one of the neighbours is debris, don't use it to calculate the
        % mean only
        cellNode_Adjacency = cellNode_Adjacency([Geo.Cells(neighboursToT1).AliveStatus] > 0);
        opposedNodes = opposedNodes([Geo.Cells(neighboursToT1).AliveStatus] > 0);
        
        % Check how they are connected (cell nodes)
        connectedNodes = cellfun(@length, cellNode_Adjacency) == length(neighboursToT1);
        nonConnectedNodes = cellfun(@length, cellNode_Adjacency) < length(neighboursToT1);
        
        % Check average distance between connected and unconnected nodes
        avgDistance_ConnNodes = mean(pdist2(Geo.Cells(ghostNodesToT1).X, vertcat(Geo.Cells(vertcat(opposedNodes{connectedNodes})).X)));
        avgDistance_NonConnNodes = mean(pdist2(Geo.Cells(ghostNodesToT1).X, vertcat(Geo.Cells(vertcat(opposedNodes{nonConnectedNodes})).X)));
        
        if avgDistance_ConnNodes < avgDistance_NonConnNodes
            continue
        end
        
        % Change the ghostPair to intercalate
        neighboursToT1 = neighboursToT1([Geo.Cells(neighboursToT1).AliveStatus] > 0);
        ghostPair = [ghostNodesToT1 neighboursToT1(connectedNodes)];
    else
        if edgeLength > Set.RemodelTol
            continue
        end
    end
    
    % Edge valence (number of shared tets)
    [valence, sharedTets] = edgeValence(Geo, ghostPair);
    
    % Face IDs involved
    faceIDs = [];
    for cellNode = sharedCellNodes'
        for numFace = 1:length(Geo.Cells(cellNode).Faces)
            face = Geo.Cells(cellNode).Faces(numFace);
            if any(ismember(ghostPair, face.ij))
                faceIDs(end+1) = face.globalIds;
            end
        end
    end
    
    % Add it to the table
    segmentFeatures(end+1, :) = table(ghostPair(1), ghostPair(2), edgeLength, valence, {sharedCellNodes}, {faceIDs}, {neighbours_1}, {neighbours_2}, {sharedTets});
end

% allTs = vertcat(Geo.Cells.T);
% cellNodes = 1:Geo.nCells;
% cellNodes = cellNodes(~cellfun(@isempty, {Geo.Cells(1:Geo.nCells).AliveStatus}));
% aliveCellNodes = cellNodes([Geo.Cells(cellNodes).AliveStatus] == 1);
% cellTets = allTs(sum(ismember(allTs, cellNodes), 2) == 1, :);
% triMesh = cellTets(sum(ismember(cellTets, aliveCellNodes), 2) == 1, :);
% triMesh = unique(sort(triMesh, 2), 'rows');
% 
% for numGhostNode = Geo.XgID
%     aspectRatio = [];
%     for tet = triMesh(any(ismember(triMesh, numGhostNode), 2), :)'
%         ghostNodes = tet(ismember(tet, Geo.XgID));
%         cellNode = tet(~ismember(tet, Geo.XgID));
%         x1 = Geo.Cells(ghostNodes(1)).X;
%         x2 = Geo.Cells(ghostNodes(2)).X;
%         x3 = Geo.Cells(ghostNodes(3)).X;
%         sideLengths(1) = norm(x1 - x2);
%         sideLengths(2) = norm(x2 - x3);
%         sideLengths(3) = norm(x3 - x1);
%         [aspectRatio(end+1)] = ComputeTriAspectRatio(sideLengths);
%     end
%     if ~isempty(aspectRatio)
%         energyPerCellAndFaces(end+1, 1:2) = horzcat(numGhostNode, mean(aspectRatio));
%     end
% end
% energyPerCellAndFaces

% for c = 1:Geo.nCells
%     for numFace = 1:length(Geo.Cells(c).Faces)
%         face = Geo.Cells(c).Faces(numFace);
%         
%         if Geo.Cells(c).AliveStatus == 1 && median([face.Tris.AspectRatio]) >= Set.RemodelTol && ~isequal(face.InterfaceType, 'CellCell') 
%             energyPerCellAndFaces(end+1, 1:6) = horzcat(c, numFace, max([face.Tris.AspectRatio]), face.globalIds, face.ij);
%         end
%     end
% end
% 

if ~isempty(segmentFeatures)
    [segmentFeatures] = sortrows(segmentFeatures, 3, 'ascend');
end

end

