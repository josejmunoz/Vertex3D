function [nodeNeighbours] = getNodeNeighbours(Geo, node, mainNode)
%GETNODENEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here

if exist('mainNode', 'var')
    allNodeTets = Geo.Cells(node).T;
    nodeNeighbours = unique(allNodeTets(any(ismember(allNodeTets, mainNode), 2), :));
    nodeNeighbours(nodeNeighbours==node) = [];
else
    nodeNeighbours = unique(Geo.Cells(node).T);
    nodeNeighbours(nodeNeighbours==node) = [];
end

end

