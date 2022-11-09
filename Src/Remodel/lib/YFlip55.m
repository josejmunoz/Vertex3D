function [Ynew, Tnew] = YFlip55(oldTets, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

Ynew = [];

Xs = unique(oldTets);
Xs_g = Xs(ismember(Xs, Geo.XgID));
Xs_c = Xs(~ismember(Xs, Geo.XgID));

% Temporary remove 4-cell tetrahedra
[tets4Cells] = get4FoldTets(Geo);
Geo = RemoveTetrahedra(Geo, tets4Cells);

[Xs_gConnectedNodes, Xs_gUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_g, Xs_g);
[Xs_cConnectedNodes, Xs_cUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_c, Xs_g);

if length(Xs_gConnectedNodes) == length(Xs_gUnconnectedNodes) && length(Xs_cConnectedNodes) == length(Xs_cUnconnectedNodes)
    if any(ismember(XsToDisconnect, Xs_gConnectedNodes)) && any(ismember(XsToDisconnect, Xs_cConnectedNodes))
        Xs_cToDisconnect = XsToDisconnect(~ismember(XsToDisconnect, Geo.XgID));
        Xs_gToDisconnect = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));
        
        XsPos(1) = Xs_gToDisconnect;
        XsPos(2) = Xs_cToDisconnect;
        XsPos(3) = setdiff(Xs_gConnectedNodes, Xs_gToDisconnect);
        XsPos(6) = setdiff(Xs_cConnectedNodes, Xs_cToDisconnect);
        XsPos(4) = Xs_gUnconnectedNodes(1);
        XsPos(5) = intersect(getNodeNeighbours(Geo, XsPos(4)), Xs_cUnconnectedNodes);
        XsPos(7) = setdiff(Xs_cUnconnectedNodes, XsPos(5));
        
        %% Connections #1: 1 mainNodes and 3 ghost node
        Tnew = [XsPos(1) XsPos(3) XsPos(4) XsPos(5)];
        
        %% Connections #2: 2 mainNodes and 2 ghost node
        Tnew = [Tnew; ...
            XsPos(7) XsPos(5) XsPos(6) XsPos(1); ...
            XsPos(7) XsPos(5) XsPos(2) XsPos(3)];
        
        %% Connections #3: 3 mainNodes and 1 ghost node
        Tnew = [Tnew; ...
            XsPos(5) XsPos(2) XsPos(4) XsPos(3); ...
            XsPos(5) XsPos(7) XsPos(1) XsPos(3)];
        
        %% Connection #4: 4 mainNodes
        Tnew(end+1, :) = [Xs_cConnectedNodes', Xs_cUnconnectedNodes'];
    else
        error('Need to check this');
    end
else
    error('Need to check this');
end

%[Ynew] = RecalculateYs(Geo, Tnew, Xs_cUnconnectedNodes, Set);
[Ynew] = RecalculateYsFromPrevious(Geo, Tnew, Xs_cUnconnectedNodes, Set);
end

