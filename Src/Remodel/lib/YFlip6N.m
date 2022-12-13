function [Ynew, Tnew] = YFlip6N(oldTets, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

Ynew = [];

ghostNodesWithoutDebris = setdiff(Geo.XgID, Geo.RemovedDebrisCells);

Xs = unique(oldTets);
Xs_g = Xs(ismember(Xs, ghostNodesWithoutDebris));
Xs_c = Xs(~ismember(Xs, ghostNodesWithoutDebris));

% Temporary remove 4-cell tetrahedra
[tets4Cells] = get4FoldTets(Geo);
Geo = RemoveTetrahedra(Geo, tets4Cells);

[Xs_gConnectedNodes, Xs_gUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_g, Xs_g(1));
[Xs_cConnectedNodes, Xs_cUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_c, Xs_g(1));

if length(Xs_gConnectedNodes) == length(Xs_gUnconnectedNodes) && length(Xs_cConnectedNodes) == length(Xs_cUnconnectedNodes)
    
    XsPos(4) = XsToDisconnect(1);
    XsPos(8) = XsToDisconnect(2);
    XsPos(3) = Xs_gUnconnectedNodes(1);
    XsPos(1) = Xs_gUnconnectedNodes(2);

    tetOf6 = oldTets(sum(ismember(oldTets, [XsPos(8) XsPos(1) XsPos(4)]), 2) > 2, :);
    tetOf6(ismember(tetOf6, [XsPos(8) XsPos(1) XsPos(4)])) = [];

    tetOf2 = oldTets(sum(ismember(oldTets, [XsPos(8) XsPos(3) XsPos(4)]), 2) > 2, :);
    tetOf2(ismember(tetOf2, [XsPos(8) XsPos(3) XsPos(4)])) = [];
    
    XsPos(2) = tetOf2;
    XsPos(6) = tetOf6;

    tetOf7 = oldTets(sum(ismember(oldTets, [XsPos(8) XsPos(2) XsPos(6)]), 2) > 2, :);
    tetOf7(ismember(tetOf7, [XsPos(8) XsPos(2) XsPos(6)])) = [];

    tetOf5 = oldTets(sum(ismember(oldTets, [XsPos(4) XsPos(2) XsPos(6)]), 2) > 2, :);
    tetOf5(ismember(tetOf5, [XsPos(4) XsPos(2) XsPos(6)])) = [];
    XsPos(7) = intersect(tetOf7, Xs_cUnconnectedNodes);
    XsPos(5) = intersect(tetOf5, Xs_cUnconnectedNodes);

    %% Connections #1: 1 mainNodes and 3 ghost node
    Tnew = [XsPos(1) XsPos(3) XsPos(8) XsPos(7); ... %3-7-8-1
        XsPos(1) XsPos(3) XsPos(4) XsPos(5)]; % 3-5-4-1

    %% Connections #2: 2 mainNodes and 2 ghost node
    Tnew = [Tnew; ...
        XsPos(7) XsPos(5) XsPos(1) XsPos(3)];

    %% Connections #3: 3 mainNodes and 1 ghost node
    Tnew = [Tnew; ...
        XsPos(7) XsPos(2) XsPos(5) XsPos(3); ...
        XsPos(5) XsPos(7) XsPos(1) XsPos(6)];

    %% Connection #4: 4 mainNodes
    Tnew(end+1, :) = [XsPos(7) XsPos(2) XsPos(5) XsPos(6)];
else
    error('Need to check this');
end
end

