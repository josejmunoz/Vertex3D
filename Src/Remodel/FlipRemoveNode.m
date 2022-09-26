function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, cellNodeLoosing, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPREMOVENODE Summary of this function goes here
%   Detailed explanation goes here

oldTets = Geo.Cells(nodeToRemove).T;
nodesToChange = getNodeNeighbours(Geo, nodeToRemove, cellNodeLoosing);
mainNodes = nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus}));
flipName = 'RemoveNode';
[Tnew, Ynew, oldTets] = ConnectTetrahedra(Geo, nodeToRemove, nodesToChange, oldTets, mainNodes, Set, flipName, cellNodeLoosing);
    
[Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);

end

