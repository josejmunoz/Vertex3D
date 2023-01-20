function [valence, sharedTets, sharedYs] = edgeValence(Geo, nodesEdge)
%EDGEVALENCE Tets in common with an edge
%   Detailed explanation goes here
nodeTets1 = sort(Geo.Cells(nodesEdge(1)).T, 2);
nodeTets2 = sort(Geo.Cells(nodesEdge(2)).T, 2);

tetIds = ismember(nodeTets1, nodeTets2, 'rows');
sharedTets = nodeTets1(tetIds, :);
sharedYs = Geo.Cells(nodesEdge(1)).Y(tetIds, :);
valence = size(sharedTets, 1);
end

