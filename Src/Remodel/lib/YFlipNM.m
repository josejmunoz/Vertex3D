function [Ynew, Tnew] = YFlipNM(oldTets, cellToIntercalateWith, oldYs, XsToDisconnect, Geo)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

allXs = vertcat(Geo.Cells.X);
Xs_gToDisconnect = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));

% Temporary remove 4-cell tetrahedra
[tets4Cells] = get4FoldTets(Geo);
Geo = RemoveTetrahedra(Geo, tets4Cells);


%% ----------- General method
visualizeTets(oldTets, vertcat(Geo.Cells.X))
%% Step 1: Keep the boundary of tets not changed.
% Or remo the edges that pass through inside the tetrahedro
boundaryNodes = unique(oldTets); % All the nodes should be boundary nodes
tris = triangulation(oldTets, allXs);

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
% possibleEdges(ismember(possibleEdges, boundaryEdges, 'rows'), :) = [];
% possibleEdges(ismember(possibleEdges, XsToDisconnect, 'rows'), :) = [];


%% Step 3: Select the edge to add
% Aim: connecting the XsToDisconnect_c to another gNode
% Based on Valence? distance? who should be intercalating with?
edgeToConnect = [cellToIntercalateWith, Xs_gToDisconnect];
possibleEdges(ismember(possibleEdges, edgeToConnect, 'rows'), :) = [];
%finalEdges = vertcat(boundaryEdges, edgeToConnect);

%% Step 4: Propagate the change to get the remaining tets
% Create tetrahedra

%% Remove edge using TetGen's algorithm
% Based on https://dl.acm.org/doi/pdf/10.1145/2629697
treeOfPossibilities = digraph;
treeOfPossibilities = addnode(treeOfPossibilities, 2);
TRemoved = {};
Tnew = {};
Ynew = {};
parentNode = 1;
arrayPos = 3;
endNode = 2;
[~, Tnew, TRemoved, treeOfPossibilities] = YFlipNM_recursive(oldTets, TRemoved, Tnew, Ynew, oldYs, Geo, possibleEdges, XsToDisconnect, treeOfPossibilities, parentNode, arrayPos);

[paths] = allpaths(treeOfPossibilities, parentNode, endNode);
newTets_tree = {};
volDiff = [];
for path =  paths'
    cPath = path{1};
    newTets = vertcat(oldTets);

    for posPath = cPath(cPath>2)
        toAdd = Tnew{posPath};
        toRemove = TRemoved {posPath};
        newTets(ismember(sort(newTets, 2), sort(toRemove, 2), 'rows'), :) = [];
        newTets = vertcat(newTets, toAdd);
    end

    itWasFound = 0;
    for numNewTet = 1:length(newTets_tree)
        newTet_tree = newTets_tree{numNewTet};
        if all(ismember(sort(newTet_tree, 2), sort(newTets, 2), 'rows'))
            itWasFound = 1;
            break
        end
    end
    if ~itWasFound
        volumes = [];
        for tet = newTets'
            [vol] = ComputeTetVolume(tet, Geo);
            volumes(end+1) = vol;
        end
        
        normVols = volumes/max(volumes);
        newTets = newTets(normVols > 0.05, :);
        newVol = sum(volumes(normVols > 0.05));

        %% Check if the volume from previous space is the same occupied by the new tets
        oldVol = 0;
        for tet = oldTets'
            [vol] = ComputeTetVolume(tet, Geo);
            oldVol = oldVol + vol;
        end

        if abs(newVol - oldVol) / oldVol <= 0.005
            newTets_tree{end+1} = newTets;
            volDiff(end+1) = abs(newVol - oldVol) / oldVol;
        end
    end
end
[~, minIndex]=min(volDiff);
Tnew = newTets_tree{minIndex};
Ynew = [];
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