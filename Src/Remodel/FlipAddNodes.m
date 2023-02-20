function [Geo, Geo_n, Geo_0, Dofs, Set, newYgIds, hasConverged] = FlipAddNodes(surroundingNodes, oldTets, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPADDNODEs Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
% Get main node (cell node)
mainNode = surroundingNodes(~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus}));
commonNodes = surroundingNodes;
% Remove the main node from the common nodes
commonNodes(ismember(commonNodes, mainNode)) = [];
flipName = 'AddNode';

if length(mainNode) > 1
    error('Too many main nodes')
end

% Add the new node in the positions (newNodes) and get the new IDs
[Geo, newNodeIDs] = AddNewNode(Geo, newNodes, commonNodes);
% Same in Geo_n
[Geo_n] = AddNewNode(Geo_n, newNodes, commonNodes);
% Same in Geo_0
[Geo_0] = AddNewNode(Geo_0, newNodes, commonNodes);

% Put together the new neighbourhood to be connected
nodesToChange = horzcat(unique(commonNodes)', newNodeIDs, mainNode);
% Connect the nodes regarding distance (delaunay method)
[Tnew] = ConnectTetrahedra(Geo, nodesToChange, oldTets, mainNode, flipName);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(tetsToChange, vertcat(Geo.Cells.X));

% Rebuild topology and run mechanics
[Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, [], oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName, [newNodeIDs -1]);
end

