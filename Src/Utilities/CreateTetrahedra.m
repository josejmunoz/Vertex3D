function [Twg] = CreateTetrahedra(trianglesConnectivity, neighboursNetwork, edgesOfVertices, xInternal, X_FaceIds, X_VerticesIds, X)
%CREATETETRAHEDRA Add connections between real nodes and ghost cells
%   Detailed explanation goes here

X_Ids = [X_FaceIds, X_VerticesIds];
Twg = [];

%% Relationships: 1 ghost node, three cell nodes
Twg_vertices = horzcat(trianglesConnectivity, X_VerticesIds');

Twg_faces = [];

Twg = vertcat(Twg_vertices, Twg_faces);

%% Relationships: 1 cell node and 3 ghost nodes
% These are the ones are with the face ghost cell on top and bottom
% 1 cell node: 1 face centre of and 2 vertices ghost nodes.
%visualizeTets(Twg(any(ismember(Twg, 1), 2), :), X)
newAdditions = [];

% Cells and faces share the same order of Ids
for numCell = xInternal'
     faceId = X_FaceIds(numCell);
     verticesToConnect = edgesOfVertices{numCell};

     newAdditions = [newAdditions; repmat([numCell, faceId], size(verticesToConnect, 1), 1), X_VerticesIds(verticesToConnect)];
end

Twg = [Twg; newAdditions];

%% Relationships: 2 ghost nodes, two cell nodes
% two of the previous ones go with 
Twg_sorted = sort(Twg(any(ismember(Twg, X_Ids), 2), :), 2);
internalNeighbourNetwork = neighboursNetwork(any(ismember(neighboursNetwork, xInternal), 2), :);
internalNeighbourNetwork = unique(sort(internalNeighbourNetwork, 2), 'rows');

newAdditions = [];
for numPair = 1:size(internalNeighbourNetwork, 1)
    [found] = ismember(Twg_sorted, internalNeighbourNetwork(numPair, :));
    newConnections = unique(Twg_sorted(sum(found, 2) == 2, 4));
    if length(newConnections) > 1
        newConnectionsPairs = nchoosek(newConnections, 2);
        newAdditions = [newAdditions; repmat(internalNeighbourNetwork(numPair, :), size(newConnectionsPairs, 1), 1), newConnectionsPairs];
    else
        error('Somewhere creating the connections and initial topology');
    end
end

Twg = [Twg; newAdditions];

end
