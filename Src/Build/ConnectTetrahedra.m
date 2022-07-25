function [newTets] = ConnectTetrahedra(Geo, newNodeIDs, oldTets)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

nodesToChange = [unique(oldTets); newNodeIDs];
newTets = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X)));

% Remove tets with all Ghost Nodes
newTets(all(ismember(newTets, Geo.XgID), 2), :) = [];
%TODO: REMOVE THE TETS THAT ADD NEW NODES TO THE CELLS
%newTets_removedNotInvolved = newTets(any(ismember(newTets, newNodeIDs), 2), :);
tetsToExclude_Possibly = newTets(~any(ismember(newTets, newNodeIDs), 2), :);
newTets_removedNotInvolved = newTets(any(ismember(newTets, newNodeIDs), 2), :);
addOrNot = [];
for newTet = tetsToExclude_Possibly'
    if mod(sum(sum(ismember(newTets_removedNotInvolved, newTet), 2) > 2), 2)
        addOrNot(end+1) = 1;
    else
        addOrNot(end+1) = 0;
    end
end

newTets = [newTets_removedNotInvolved; tetsToExclude_Possibly(addOrNot==1, :)];
end

