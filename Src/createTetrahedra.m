function [Twg] = createTetrahedra(trianglesConnectivity, neighboursNetwork, verticesNodesNetwork, X, xInternal, X_bottomNodes, X_bottomFaceIds, X_bottomVerticesIds)
%CREATETETRAHEDRA Add connections between real nodes and ghost cells
%   Detailed explanation goes here

X_bottomIds = [X_bottomFaceIds, X_bottomVerticesIds];

% Relationships: 1 ghost node, three cell nodes
verticesTriangles_1 = X(trianglesConnectivity(:, 1), :);
verticesTriangles_2 = X(trianglesConnectivity(:, 2), :);
verticesTriangles_3 = X(trianglesConnectivity(:, 3), :);
meanTriangles = arrayfun(@(x, y, z) mean([x, y, z]), verticesTriangles_1, verticesTriangles_2, verticesTriangles_3);

[~, indices] = pdist2(X_bottomNodes, meanTriangles, 'euclidean', 'smallest', 1);
Twg = horzcat(trianglesConnectivity, X_bottomIds(indices)');

% Relationships: 2 ghost nodes, two cell nodes
% two of the previous ones go with 
Twg_sorted = sort(Twg(any(ismember(Twg, X_bottomIds), 2), :), 2);
internalNeighbourNetwork = neighboursNetwork(any(ismember(neighboursNetwork, xInternal), 2), :);
internalNeighbourNetwork = unique(sort(internalNeighbourNetwork, 2), 'rows');

newAdditions = [];
for numPair = 1:size(internalNeighbourNetwork, 1)
    [found] = ismember(Twg_sorted, internalNeighbourNetwork(numPair, :));
    newConnections = unique(Twg_sorted(sum(found, 2) == 2, 4));
    if length(newConnections) == 2
        newAdditions = [newAdditions; internalNeighbourNetwork(numPair, :) newConnections'];
    end
end

Twg = [Twg; newAdditions];

% Relationships: 1 cell node and 3 ghost nodes
% These are the ones are with the face ghost cell on top and bottom
% 1 cell node: 1 face centre of and 2 vertices ghost nodes.

newAdditions = [];

% Cells and faces share the same order of Ids
for numCell = xInternal
     faceId = X_bottomFaceIds(numCell);
     verticesToConnect = verticesNodesNetwork(verticesNodesNetwork(:, 1) == numCell, 2);
     
     newAdditions = [newAdditions; [numCell, faceId]];
end


end

