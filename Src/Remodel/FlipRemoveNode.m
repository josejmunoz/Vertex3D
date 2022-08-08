function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPREMOVENODE Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;

oldTets = Geo.Cells(nodeToRemove).T;
nodesToChange = getNodeNeighbours(Geo, nodeToRemove);
mainNode = nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus}));

if length(mainNode) > 1
    return
end

[Tnew] = ConnectTetrahedra(Geo, nodesToChange, oldTets, mainNode);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));

[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'RemoveNode');

end

