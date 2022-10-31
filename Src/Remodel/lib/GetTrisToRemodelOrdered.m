function [segmentFeatures] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here
possibleGhostPairs = nchoosek(Geo.XgID, 2);

segmentFeatures = table();
for ghostPair = possibleGhostPairs'
    if ismember(ghostPair(2), Geo.Cells(ghostPair(1)).T)
        % Edge length
        x1 = Geo.Cells(ghostPair(1)).X;
        x2 = Geo.Cells(ghostPair(2)).X;
        edgeLength = norm(x2 - x1);
        
        % Edge valence (number of shared tets)
        [valence, sharedTets] = edgeValence(Geo, ghostPair);
        
        % Number of cell nodes shared
        sharedCellNodes = unique(sharedTets(~ismember(sharedTets, Geo.XgID)));
        
        segmentFeatures(end+1, :) = table(ghostPair(1), ghostPair(2), edgeLength, valence, {sharedCellNodes});
    end
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
    [segmentFeatures] = sortrows(segmentFeatures, 2, 'descend');
end

end

