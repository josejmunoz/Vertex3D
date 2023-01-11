function [Ynew, Tnew] = YFlip55(oldTets, XsToDisconnect, Geo, Set)
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


if length(Xs_cConnectedNodes) == length(Xs_c) && length(Xs_gConnectedNodes) == length(Xs_gUnconnectedNodes)
    Xs_cToDisconnect = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));
    Xs_gToDisconnect = XsToDisconnect(~ismember(XsToDisconnect, Geo.XgID));
    
    [Xs_cConnectedNodes, Xs_cUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_g, Xs_g(1));
    [Xs_gConnectedNodes, Xs_gUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_c, Xs_g(1));
    
    
    Xs_g = Xs(~ismember(Xs, ghostNodesWithoutDebris));
    Xs_c = Xs(ismember(Xs, ghostNodesWithoutDebris));
else
    Xs_cToDisconnect = XsToDisconnect(~ismember(XsToDisconnect, Geo.XgID));
    Xs_gToDisconnect = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));
end 

if length(Xs_gConnectedNodes) == length(Xs_g) && length(Xs_cConnectedNodes) == length(Xs_cUnconnectedNodes)
    if any(ismember(XsToDisconnect, Xs_gConnectedNodes)) && any(ismember(XsToDisconnect, Xs_cConnectedNodes))
        XsPos(1) = Xs_gToDisconnect;
        XsPos(2) = Xs_cToDisconnect;
        
        Xs_gRemaining = setdiff(Xs_gConnectedNodes, Xs_gToDisconnect);
        
        %% COULDNT FIND A BETTER WAY OF GETTING WHICH NODE SHOULD BE THE 3 AND 4
        aliveStatus = [Geo.Cells(Xs_gRemaining).AliveStatus] == 1;
        if any(~aliveStatus) && ~all(aliveStatus == 0)
            XsPos(3) = Xs_gRemaining(aliveStatus);
            XsPos(4) = Xs_gRemaining(~aliveStatus);
        else
            if norm(Geo.Cells(Xs_gRemaining(2)).X(1:2) - Geo.Cells(XsPos(2)).X(1:2)) > norm(Geo.Cells(Xs_gRemaining(1)).X(1:2) - Geo.Cells(XsPos(2)).X(1:2))
                XsPos(3) = Xs_gRemaining(1);
                XsPos(4) = Xs_gRemaining(2);
            else
                XsPos(3) = Xs_gRemaining(2);
                XsPos(4) = Xs_gRemaining(1);
            end
        end
        XsPos(6) = setdiff(Xs_cConnectedNodes, Xs_cToDisconnect);
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
        %Tnew(end+1, :) = [Xs_cConnectedNodes', Xs_cUnconnectedNodes'];
    else
        error('Need to check this');
    end
else
    error('Need to check this');
end
end

