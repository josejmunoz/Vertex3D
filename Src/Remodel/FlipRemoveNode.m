function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, internalFlip)
%FLIPREMOVENODE Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;

oldTets = Geo.Cells(nodeToRemove).T;
nodesToChange = getNodeNeighbours(Geo, nodeToRemove);
mainNodes = nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus}));

Tnew = [];
for mainNode = mainNodes'
    currentOldTets = oldTets(any(ismember(oldTets, mainNode), 2), :);
    currenteNodesToChange = unique(currentOldTets);
    currenteNodesToChange(nodeToRemove == currenteNodesToChange) = [];
    [c_Tnew] = ConnectTetrahedra(Geo, currenteNodesToChange, currentOldTets, mainNode);
    Tnew = vertcat(Tnew, c_Tnew);
end

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));

if length(mainNodes) > 1
    internalFlip = 1;
else
    internalFlip = 0;
end
    
[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'RemoveNode', internalFlip);

end

