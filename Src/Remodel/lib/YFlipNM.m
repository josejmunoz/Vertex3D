function [Ynew, Tnew] = YFlipNM(oldTets, cellToIntercalateWith, oldYs, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

oldTets_original = oldTets;
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
tris = triangulation(oldTets, allXs);
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
TRemoved = {[]};
Tnew = {[]};
Ynew = {[]};
treeDepth = 1;
[Ynew, Tnew, TRemoved] = YFlipNM_recursive(oldTets, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect, treeDepth);
Tnew
TRemoved

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