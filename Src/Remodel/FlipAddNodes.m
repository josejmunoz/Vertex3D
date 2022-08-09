function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNodes(surroundingNodes, tetsToChange, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPADDNODEs Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
mainNode = surroundingNodes(~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus}));
commonNodes = surroundingNodes;
commonNodes(ismember(commonNodes, mainNode)) = [];

if length(mainNode) > 1
    return
end

[Geo, newNodeIDs] = AddNewNode(Geo, newNodes);
[Geo_n] = AddNewNode(Geo_n, newNodes);

nodesToChange = horzcat(unique(commonNodes)', newNodeIDs, mainNode);
[Tnew] = ConnectTetrahedra(Geo, nodesToChange, tetsToChange, mainNode);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));

[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, tetsToChange, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'AddNode');
end

