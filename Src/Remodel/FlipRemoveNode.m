function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, cellNodeLoosing, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPREMOVENODE Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;

oldTets = Geo.Cells(nodeToRemove).T;
nodesToChange = getNodeNeighbours(Geo, nodeToRemove, cellNodeLoosing);
mainNodes = nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus}));
flipName = 'RemoveNode';
[Tnew, Ynew] = ConnectTetrahedra(Geo, nodeToRemove, nodesToChange, oldTets, mainNodes, Set, flipName, cellNodeLoosing);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));

% if length(mainNodes) > 1
%     flipName = 'RemoveNode: 2+ Main nodes';
% else
%     flipName = 'RemoveNode';
% end
    
[Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);

end

