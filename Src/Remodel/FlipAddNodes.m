function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNodes(surroundingNodes, tetsToChange, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPADDNODEs Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
% Get main node (cell node)
mainNode = surroundingNodes(~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus}));
commonNodes = surroundingNodes;
% Remove the main node from the common nodes
commonNodes(ismember(commonNodes, mainNode)) = [];

if length(mainNode) > 1
    return
end

% Add the new node in the positions (newNodes) and get the new IDs
[Geo, newNodeIDs] = AddNewNode(Geo, newNodes);
% Same in Geo_n
[Geo_n] = AddNewNode(Geo_n, newNodes);

% Put together the new neighbourhood to be connected
nodesToChange = horzcat(unique(commonNodes)', newNodeIDs, mainNode);
% Connect the nodes regarding distance (delaunay method)
[Tnew] = ConnectTetrahedra(Geo, nodesToChange, tetsToChange, mainNode);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));

% Rebuild topology and run mechanics
[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, tetsToChange, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'AddNode');
end

