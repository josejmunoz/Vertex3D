function [Ynew, Tnew] = YFlipNM(oldTets, cellToIntercalateWith, oldYs, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

Ynew = [];
allXs = vertcat(Geo.Cells.X);
ghostNodesWithoutDebris = setdiff(Geo.XgID, Geo.RemovedDebrisCells);

Xs = unique(oldTets);
Xs_g = Xs(ismember(Xs, ghostNodesWithoutDebris));
Xs_c = Xs(~ismember(Xs, ghostNodesWithoutDebris));
Xs_cToDisconnect = XsToDisconnect(~ismember(XsToDisconnect, Geo.XgID));
Xs_gToDisconnect = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));
intercalationFlip = 0;
if length(Xs_c) == 4
    intercalationFlip = 1;
end

% Temporary remove 4-cell tetrahedra
[tets4Cells] = get4FoldTets(Geo);
Geo = RemoveTetrahedra(Geo, tets4Cells);

[Xs_gConnectedNodes, Xs_gUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_g, Xs_g(1));
[Xs_cConnectedNodes, Xs_cUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_c, Xs_g(1));

%% ----------- General method
visualizeTets(oldTets, vertcat(Geo.Cells.X))
%% Step 1: Keep the boundary of tets not changed.
% Or remo the edges that pass through inside the tetrahedro
boundaryNodes = unique(oldTets); % All the nodes should be boundary nodes
tris = triangulation(oldTets, vertcat(Geo.Cells.X));
trisEdges = tris.edges;

% nodes should be connected with two of the same layer.
tets_Ghost = tris.ConnectivityList(sum(ismember(tris.ConnectivityList, Geo.XgID), 2)>2, :);
tets_Ghost = sort(tets_Ghost, 2);
tris_Ghost = tets_Ghost(:, 2:end);
if any(~ismember(tris_Ghost(:), Geo.XgID))
    error('Yflip55')
end
tets_Cells = tris.ConnectivityList(sum(~ismember(tris.ConnectivityList, Geo.XgID), 2)>2, :);
tets_Cells = sort(tets_Cells, 2);
tris_Cells = tets_Cells(:, 1:3);
if any(ismember(tris_Cells(:), Geo.XgID))
    error('Yflip55')
end
[nodesExt_Cells, boundary_Ghost] = GetBoundary2D(tris_Ghost, allXs);
[nodesExt_Ghost, boundary_Cells] = GetBoundary2D(tris_Cells, allXs);
boundaryEdges = vertcat(sort(boundary_Ghost, 2), sort(boundary_Cells, 2));

%% Step 2: Get edges that can be added when removing the other one
% Find the edges that may be connected and are not connected
possibleEdges = nchoosek(boundaryNodes, 2);
possibleEdges(ismember(possibleEdges, boundaryEdges, 'rows'), :) = [];
possibleEdges(ismember(possibleEdges, XsToDisconnect, 'rows'), :) = [];

%Remove other impossible edges


%% Step 3: Select the edge to add
% Aim: connecting the XsToDisconnect_c to another gNode
% Based on Valence? distance? who should be intercalating with?
edgeToConnect = [cellToIntercalateWith, Xs_gToDisconnect];
possibleEdges(ismember(possibleEdges, edgeToConnect, 'rows'), :) = [];
finalEdges = vertcat(boundaryEdges, edgeToConnect);

possibleEdges

%% Step 4: Propagate the change to get the remaining tets
% Create tetrahedra

%% Remove edge using TetGen's algorithm
% Based on https://dl.acm.org/doi/pdf/10.1145/2629697
if size(oldTets, 1) == 3

else
    % https://link.springer.com/article/10.1007/s00366-016-0480-z#Fig2
    possibleEdgesToKeep = [];
    for numPair = 1:size(possibleEdges, 1)
        [valence, sharedTets, tetIds] = edgeValenceT(oldTets, possibleEdges(numPair, :));
        
        % Valence == 1, is an edge that can be removed.
        % Valence == 2, a face can be removed.
        if valence == 2
            possibleEdgesToKeep(end+1, :) = possibleEdges(numPair, :);
            [Ynew, Tnew_23] = YFlip23(oldYs, oldTets, tetIds, Geo);
            Tnew_23
            [Ynew, Tnew_32] = YFlip32(oldYs, Tnew_23, [1 2 3], Geo);
            %% NOW WHAT???
            % What I've understand is, that you should do a series of
            % flip2-3 and a flip 3-2 to remove the edge at the end.
        end
    end
end

%% By chatGPT
% Flip edges adjacent to the target edge until it is no longer adjacent to any tetrahedra
while ~isempty(adjacent_tetrahedra)
    % Choose a random adjacent tetrahedron
    tet = adjacent_tetrahedra(1);
    
    % Find the other vertex of the tetrahedron
    other_vertex = setdiff(T(tet,:), edge);
    
    % Find the edges opposite to the other vertex
    opposite_edges = nchoosek(T(tet,:), 3);
    opposite_edges = opposite_edges(sum(ismember(opposite_edges, other_vertex), 2) == 2, :);
    
    % Flip each opposite edge that is adjacent to another tetrahedron
    for i = 1:size(opposite_edges, 1)
        opposite_vertex = setdiff(opposite_edges(i,:), [other_vertex, edge]);
        adjacent_tetrahedra = setdiff(adjacent_tetrahedra, find(sum(ismember(T, opposite_edges(i,:)), 2) == 2));
        if ~ismember(opposite_vertex, edge)
            [V, T] = flip_nm(V, T, tet, find(sum(ismember(T, opposite_edges(i,:)), 2) == 2));
        end
    end
    
    % Update the set of adjacent tetrahedra
    adjacent_tetrahedra = setdiff(adjacent_tetrahedra, tet);
end

% Remove the tetrahedra that were adjacent to the target edge
T_new = T(~sum(ismember(T, edge), 2), :);

% Update the vertex coordinates
V_new = V;

end

function [nodesExt, pairsExt]=GetBoundary2D(T,X)
np=size(X,1);
nele=size(T,1);
nodesExt=zeros(1,np);
pairsExt=[];
for e=1:nele
    Te=[T(e,:) T(e,1)];
    Sides=[0 0 0];
    for s=1:3
        n=Te(s:s+1);
        for d=1:nele
            if sum(ismember(n,T(d,:)))==2 && d~=e
                Sides(s)=1;
                break;
            end
        end
        if Sides(s)==0
            nodesExt(Te(s:s+1))=Te(s:s+1);
            pairsExt(end+1, 1:2) = Te(s:s+1);
        end
    end
end
nodesExt(nodesExt==0)=[];
end