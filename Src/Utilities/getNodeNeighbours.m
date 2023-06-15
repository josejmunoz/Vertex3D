function [nodeNeighbours] = getNodeNeighbours(Geo, node, mainNode)
%GETNODENEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here

if exist('mainNode', 'var')
    allNodeTets = vertcat(Geo.Cells(node).T);
    nodeNeighbours = unique(allNodeTets(any(ismember(allNodeTets, mainNode), 2), :));
else
    nodeNeighbours = unique(vertcat(Geo.Cells(node).T));
end

nodeNeighbours(ismember(nodeNeighbours, node)) = [];

end