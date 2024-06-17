function [Ynew, Tnew] = YFlipNM(oldTets, cellToIntercalateWith, oldYs, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

Xs_gToDisconnect = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));

% Temporary remove 4-cell tetrahedra
[tets4Cells] = get4FoldTets(Geo);
Geo = RemoveTetrahedra(Geo, tets4Cells);
tets4Cells = unique(sort(tets4Cells, 2), 'rows');
ghostNodesWithoutDebris = setdiff(Geo.XgID, Geo.RemovedDebrisCells);

Xs = unique(oldTets);
Xs_c = Xs(~ismember(Xs, ghostNodesWithoutDebris));
intercalationFlip = 0;
if length(Xs_c) == 4
    intercalationFlip = 1;
end

%% ----------- General method
%% Step 1: Keep the boundary of tets not changed.
% Or remo the edges that pass through inside the tetrahedro
boundaryNodes = unique(oldTets); % All the nodes should be boundary nodes

%% Step 2: Get edges that can be added when removing the other one
% Find the edges that may be connected and are not connected
possibleEdges = nchoosek(boundaryNodes, 2);

%% Step 3: Select the edge to add
% Aim: connecting the XsToDisconnect_c to another gNode
% Based on Valence? distance? who should be intercalating with?
edgeToConnect = sort([cellToIntercalateWith, Xs_gToDisconnect]);
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
vol_diff = Inf;
cell_winning = -Inf;
new_tets = [];
for path =  paths'
    cPath = path{1};
    newTets = vertcat(oldTets);

    for posPath = cPath(cPath>2)
        toAdd = Tnew{posPath};
        toRemove = TRemoved {posPath};
        newTets(ismember(sort(newTets, 2), sort(toRemove, 2), 'rows'), :) = [];
        newTets = vertcat(newTets, toAdd);
    end

    %No tets with 4 ghost nodes
    %if all(sum(ismember(newTets, Geo.XgID), 2) < 4)
    current_won_valence = sum(any(ismember(newTets, cellToIntercalateWith), 2))/size(newTets, 1);
    if current_won_valence >= cell_winning
        volumes = [];
        for tet = newTets'
            [vol] = ComputeTetVolume(tet, Geo);
            volumes(end+1) = vol;
        end
        
        newVol = sum(volumes);

        %% Check if the volume from previous space is the same occupied by the new tets
        oldVol = 0;
        for tet = oldTets'
            [vol] = ComputeTetVolume(tet, Geo);
            oldVol = oldVol + vol;
        end
        
        current_vol_diff = abs(newVol - oldVol) / oldVol;

        if current_vol_diff < vol_diff
            try
                if intercalationFlip
                    Xs_c = Xs(~ismember(Xs, ghostNodesWithoutDebris));
                    if ~ismember(sort(Xs_c)', sort(newTets, 2), 'rows')
                        newTets(end+1, :) = Xs_c;
                    end
                end
                [Geo_new] = RemoveTetrahedra(Geo, oldTets);
                Set.TryingFlips = 1;
                [Geo_new] = AddTetrahedra(Geo_new, Geo, [newTets; tets4Cells], [], Set);
                Set.TryingFlips = 0;
                Rebuild(Geo_new, Set);
                new_tets = newTets;
                cell_winning = current_won_valence
                vol_diff = current_vol_diff
            catch ex
                disp(newTets)
                disp(ex.message)
            end
        end
    end
    %end
end
Tnew = new_tets;
Ynew = [];
end