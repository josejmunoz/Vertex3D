function [Ynew, Tnew] = YFlip7N(oldTets, XsToDisconnect, Geo, Set)
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

if length(Xs_gConnectedNodes) == 1 && length(Xs_cConnectedNodes) == length(Xs_cUnconnectedNodes)
    if any(ismember(XsToDisconnect, Xs_gConnectedNodes)) && any(ismember(XsToDisconnect, Xs_cConnectedNodes))
        Xs_cToDisconnect = XsToDisconnect(~ismember(XsToDisconnect, ghostNodesWithoutDebris));
        Xs_gToDisconnect = XsToDisconnect(ismember(XsToDisconnect, ghostNodesWithoutDebris));
        
        XsPos(2) = Xs_cToDisconnect;
        XsPos(1) = Xs_gToDisconnect;
        
        XsPos(6) = setdiff(Xs_cConnectedNodes, XsPos(2));
        
        remainingCells = setdiff(Xs_c, XsPos);
        remainingGhost = setdiff(Xs_g, Xs_gToDisconnect);
        
        XsPos(5) = remainingCells(1);
        XsPos(7) = remainingCells(2);
        
        XsPos(4) = intersect(oldTets(any(ismember(oldTets, XsPos(5)), 2), :), remainingGhost);
        XsPos(9) = intersect(oldTets(any(ismember(oldTets, XsPos(7)), 2), :), remainingGhost);
        
        XsPos(3) = setdiff(oldTets(any(ismember(oldTets, XsPos(4)), 2), :), XsPos);
        XsPos(8) = setdiff(oldTets(any(ismember(oldTets, XsPos(9)), 2), :), XsPos);
        
        %% Option 1
        Tnew = [XsPos(3) XsPos(4) XsPos(8) XsPos(2); ...
            XsPos(9) XsPos(8) XsPos(4) XsPos(7); ...
            XsPos(5) XsPos(6) XsPos(7) XsPos(1); ...
            XsPos(9) XsPos(1) XsPos(4) XsPos(7)]
        
        %% Option 2
        Tnew = [XsPos(1) XsPos(3) XsPos(4) XsPos(5); ...
            XsPos(3) XsPos(8) XsPos(1) XsPos(7); ...
            XsPos(9) XsPos(8) XsPos(1) XsPos(7); ...
            XsPos(7) XsPos(6) XsPos(5) XsPos(1); ...
            XsPos(7) XsPos(2) XsPos(5) XsPos(3); ...
            XsPos(7) XsPos(1) XsPos(5) XsPos(3)];
        
        Tnew(end+1, :) = [XsPos(2) XsPos(5) XsPos(7) XsPos(6)];
        
    else
        disp('Need to check this');
        Tnew = [];
    end
elseif length(Xs_cConnectedNodes) == length(Xs_c) 
    disp('Need to check this');
    Tnew = [];
else
    disp('Need to check this');
    Tnew = [];
end
end