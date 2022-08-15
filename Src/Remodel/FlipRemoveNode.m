function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPREMOVENODE Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;

oldTets = Geo.Cells(nodeToRemove).T;
nodesToChange = getNodeNeighbours(Geo, nodeToRemove);
mainNodes = nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus}));

if length(mainNodes) > 1
    flipName = 'RemoveNode: 2+ Main nodes';
else
    flipName = 'RemoveNode';
end

Tnew = [];
for mainNode = mainNodes'
    currentOldTets = oldTets(any(ismember(oldTets, mainNode), 2), :);
    currenteNodesToChange = unique(currentOldTets);
    currenteNodesToChange(nodeToRemove == currenteNodesToChange) = [];
    [c_Tnew] = ConnectTetrahedra(Geo, currenteNodesToChange, currentOldTets, mainNode, flipName);
    Tnew = vertcat(Tnew, c_Tnew);
    if CheckOverlappingTets(currentOldTets, c_Tnew, Geo, flipName)
        Tnew = [];
        return
    end
end

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));
    
[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);

end

