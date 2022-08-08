function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNode(surroundingNodes, tetsToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPADDNODE Summary of this function goes here
%   Detailed explanation goes here


mainNode = surroundingNodes(~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus}));
commonNodes = surroundingNodes;
commonNodes(ismember(commonNodes, mainNode)) = [];

if length(mainNode) > 1
    return
end

[Geo, newNodeIDs] = AddNewNode(Geo, mean(vertcat(Geo.Cells(commonNodes).X)));
[Geo_n] = AddNewNode(Geo_n, mean(vertcat(Geo.Cells(commonNodes).X)));

nodesToChange = [unique(commonNodes)'; newNodeIDs];
[Tnew] = ConnectTetrahedra(Geo, nodesToChange, tetsToChange, mainNode);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));

[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, tetsToChange, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'AddNode');
end

