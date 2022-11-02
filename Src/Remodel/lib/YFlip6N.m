function [Ynew, Tnew] = YFlip6N(oldYs, oldTets, XsToDisconnect, Geo, Set)
%YFLIP6N Summary of this function goes here
%   Detailed explanation goes here
Tnew = [];
Ynew = [];

Xs = unique(oldTets);
Xs_g = Xs(ismember(Xs, Geo.XgID));
Xs_c = Xs(~ismember(Xs, Geo.XgID));

[Xs_gConnectedNodes, Xs_gUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_g);
[Xs_cConnectedNodes, Xs_cUnconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs_c);

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
        XsPos(8) = Xs_gUnconnectedNodes(2);
        XsPos(7) = intersect(getNodeNeighbours(Geo, XsPos(8)), Xs_cUnconnectedNodes);
        
        %% Connections #1: 1 mainNodes and 3 ghost node
        Tnew = [XsPos(1) XsPos(3) XsPos(8) XsPos(7); ... %3-7-8-1
            XsPos(1) XsPos(3) XsPos(4) XsPos(5)]; % 3-5-4-1
        
        %% Connections #2: 2 mainNodes and 2 ghost node
        Tnew = [Tnew; ...
            XsPos(7) XsPos(5) XsPos(6) XsPos(1); ...
            XsPos(7) XsPos(5) XsPos(2) XsPos(3)];
        
        %% Connections #3: 3 mainNodes and 1 ghost node
        Tnew = [Tnew; ...
            XsPos(7) XsPos(2) XsPos(8) XsPos(3); ...
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

mainNodesToConnect = Xs_cUnconnectedNodes;
nodesChanged = unique(Tnew(:));
oldTets = oldTets(sum(ismember(oldTets, nodesChanged), 2) > 3, :);

%% Recalculate Ys
allTs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).T);
allYs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).Y);
for numTet = 1:size(Tnew, 1)
    tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 2;
    
    mainNode_current = mainNodesToConnect(ismember(mainNodesToConnect, Tnew(numTet, :)));
    if any(tetsToUse)
        contributionOldYs = 1;
        Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * ComputeY(vertcat(Geo.Cells(Tnew(numTet, :)).X), Geo.Cells(mainNode_current(1)).X, length([Geo.Cells(Tnew(numTet, :)).AliveStatus]), Set);
    else
        contributionOldYs = 0.9;
        tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 1;
        
        if any(ismember(Tnew(numTet, :), Geo.XgTop))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
        elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
            tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
        end
        
        if any(tetsToUse)
            Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * ComputeY(vertcat(Geo.Cells(Tnew(numTet, :)).X), Geo.Cells(mainNode_current(1)).X, length([Geo.Cells(Tnew(numTet, :)).AliveStatus]), Set);
        else
            error('Need to check this!');
        end
    end
            end


end

