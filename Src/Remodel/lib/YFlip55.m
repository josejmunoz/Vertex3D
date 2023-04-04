function [Ynew, Tnew] = YFlip55(oldTets, oldYs, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here

Ynew = [];

ghostNodesWithoutDebris = setdiff(Geo.XgID, Geo.RemovedDebrisCells);

Xs = unique(oldTets);
Xs_g = Xs(ismember(Xs, ghostNodesWithoutDebris));
Xs_c = Xs(~ismember(Xs, ghostNodesWithoutDebris));
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
tris.edges

% nodes should be connected with two of the same layer and 1 different in
% the boundary.

boundaryEdges = [];

%% Step 2: Get edges that can be added when removing the other one
% Find the edges that may be connected and are not connected
possibleEdges = nchoosek(boundaryNodes, 2);
possibleEdges = setdiff(possibleEdges, boundaryEdges);
possibleEdges = setdiff(possibleEdges, XsToDisconnect);

%Remove other impossible edges

%% Step 3: Select the edge to add
% Based on connecting the 

%% Step 4: Propagate the change to get the remaining tets


%% ----------- Specific method
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
        if all(ismember(Xs_gRemaining, Geo.XgID) == 0)
            aliveStatus = [Geo.Cells(Xs_gRemaining).AliveStatus] == 1;
        else
            aliveStatus = 0;
        end
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
        if intercalationFlip
            Xs_c = Xs(~ismember(Xs, ghostNodesWithoutDebris));
            Tnew(end+1, :) = Xs_c;
        end
        %% Recalculate Ys: the new way.
        Tnew_cNodes = sum(~ismember(Tnew, Geo.XgID), 2);
        oldTets_cNodes = sum(~ismember(oldTets, Geo.XgID), 2);
%         Ynew = RecalculateYsFromPrevious(Geo, Tnew, Xs, Set);
%         for nTet = 1:size(Tnew, 1)
%             cTet = Tnew(nTet, :);
%             Ynew(nTet, :) = mean(oldYs(sum(ismember(oldTets, cTet), 2) > 2 & Tnew_cNodes(nTet) == oldTets_cNodes, :), 1);
%         end
    else
        %error('Need to check this');
        Tnew = [];
    end
else
    %error('Need to check this');
    Tnew = [];
end
end

function nodesExt=GetBoundary3D(T,X)
np=size(X,1);
nele=size(T,1);
nodesExt=zeros(1,np);
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
        end
    end
end
nodesExt(nodesExt==0)=[];
end