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