function [nodeNeighbours] = getNodeNeighbours(Geo, node)
%GETNODENEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here

nodeNeighbours = unique(Geo.Cells(node).T);
nodeNeighbours(nodeNeighbours==node) = [];

end

