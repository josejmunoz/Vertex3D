function [Twg] = createTetrahedra(trianglesConnectivity, neighboursNetwork, edgesOfVertices, X, xInternal, X_Nodes, X_FaceIds, X_VerticesIds)
%CREATETETRAHEDRA Add connections between real nodes and ghost cells
%   Detailed explanation goes here

X_Ids = [X_FaceIds, X_VerticesIds];

% Relationships: 1 ghost node, three cell nodes
verticesTriangles_1 = X(trianglesConnectivity(:, 1), :);
verticesTriangles_2 = X(trianglesConnectivity(:, 2), :);
verticesTriangles_3 = X(trianglesConnectivity(:, 3), :);
meanTriangles = arrayfun(@(x, y, z) mean([x, y, z]), verticesTriangles_1, verticesTriangles_2, verticesTriangles_3);

[~, indices] = pdist2(X_Nodes, meanTriangles, 'euclidean', 'smallest', 1);
Twg = horzcat(trianglesConnectivity, X_Ids(indices)');

% Relationships: 2 ghost nodes, two cell nodes
% two of the previous ones go with 
Twg_sorted = sort(Twg(any(ismember(Twg, X_Ids), 2), :), 2);
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
for numCell = xInternal'
     faceId = X_FaceIds(numCell);
     verticesToConnect = edgesOfVertices{numCell};
     verticesToConnect = verticesToConnect(1:end-1, :);
     
     newAdditions = [newAdditions; repmat([numCell, faceId], size(verticesToConnect, 1), 1), X_VerticesIds(verticesToConnect)];
end

Twg = [Twg; newAdditions];
end

