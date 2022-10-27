function [valence, sharedTets] = edgeValence(Geo, nodesEdge)
%EDGEVALENCE Tets in common with an edge
%   Detailed explanation goes here
nodeTets1 = sort(Geo.Cells(nodesEdge(1)).T, 2);
nodeTets2 = sort(Geo.Cells(nodesEdge(2)).T, 2);

sharedTets = nodeTets1(ismember(nodeTets1, nodeTets2, 'rows'), :);
valence = size(sharedTets, 1);
end

