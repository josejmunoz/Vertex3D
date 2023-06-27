function [Ynew, Tnew] = YFlip44(oldTets, XsToDisconnect, Geo, Set)
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

if any(ismember(XsToDisconnect, Xs_gConnectedNodes)) && any(ismember(XsToDisconnect, Xs_cConnectedNodes))
    nodeToDisconnect_c = XsToDisconnect(~ismember(XsToDisconnect, Geo.XgID));
    nodeToDisconnect_g = XsToDisconnect(ismember(XsToDisconnect, Geo.XgID));
    if length(Xs_gConnectedNodes) == length(Xs_gUnconnectedNodes) && length(Xs_cConnectedNodes) == length(Xs_c)
        %% XsToDisconnect cell should be connected to all of them
        % Connect the other cell to all of them
        nodeToConnect_c = setdiff(Xs_cConnectedNodes, XsToDisconnect);
        nodeCellNeighbours = getNodeNeighboursPerDomain(Geo, nodeToConnect_c, Xs_g(1));
        nodeCellNeighbours_g = nodeCellNeighbours(ismember(nodeCellNeighbours, Geo.XgID));
        nodeToConnect_g = setdiff(Xs_g, nodeCellNeighbours_g);
        nodeToConnect_g = intersect(Xs_gUnconnectedNodes, nodeToConnect_g);
        
        if length(nodeToConnect_g) > 1
            disp('Probleeeemmmmm')
            Tnew = [];
            return
        end
        
        XsPos = [-1 -1 nodeToConnect_g nodeToDisconnect_g nodeToConnect_c nodeToDisconnect_c];
        XsPos(1) = setdiff(Xs_gUnconnectedNodes, nodeToConnect_g);
        XsPos(2) = setdiff(Xs_gConnectedNodes, nodeToDisconnect_g);
        
        Tnew = [XsPos(1) XsPos(3) XsPos(4) XsPos(5); ...
            XsPos(1) XsPos(3) XsPos(2) XsPos(5); ...
            XsPos(1) XsPos(2) XsPos(5) XsPos(6)]; %...
           % XsPos(1) XsPos(5) XsPos(4) XsPos(6)];
    elseif length(Xs_gConnectedNodes) == length(Xs_g) && length(Xs_cConnectedNodes) == length(Xs_c)
        nodeToConnect_c = setdiff(Xs_cConnectedNodes, XsToDisconnect);
        nodeToConnect_c = nodeToConnect_c([Geo.Cells(nodeToConnect_c).AliveStatus] == 1);
        
        if length(nodeToConnect_c) > 1
            disp('Probleeeemmmmm')
            Tnew = [];
            return
        end
        
        nodeCellNeighbours = getNodeNeighboursPerDomain(Geo, nodeToConnect_c, Xs_g(1));
        nodeCellNeighbours_g = nodeCellNeighbours(ismember(nodeCellNeighbours, Geo.XgID));
        nodeToConnect_g = setdiff(Xs_g, nodeCellNeighbours_g);
        nodeToConnect_g = intersect(Xs_gConnectedNodes, nodeToConnect_g);
        
        XsPos = [-1 -1 nodeToConnect_g nodeToDisconnect_g nodeToConnect_c nodeToDisconnect_c];
        XsPos(1) = setdiff(Xs_g, [nodeToConnect_g nodeToDisconnect_g]);
        XsPos(2) = setdiff(Xs_c, [nodeToConnect_c nodeToDisconnect_c]);
        
        Tnew = [XsPos(2) XsPos(5) XsPos(6) XsPos(3); ...
            XsPos(5) XsPos(1) XsPos(3) XsPos(4); ...
            XsPos(5) XsPos(6) XsPos(3) XsPos(1); ...
            XsPos(5) XsPos(2) XsPos(3) XsPos(4)];
    else
        Tnew = [];
    end
else
    Tnew = [];
end
end