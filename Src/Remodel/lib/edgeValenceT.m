function [valence, sharedTets, tetIds] = edgeValenceT(tets, nodesEdge)
%EDGEVALENCE Tets in common with an edge
%   Detailed explanation goes here
tets1 = tets(any(ismember(tets, nodesEdge(1)), 2), :);
tets2 = tets(any(ismember(tets, nodesEdge(2)), 2), :);

nodeTets1 = sort(tets1, 2);
nodeTets2 = sort(tets2, 2);

tetIds = ismember(nodeTets1, nodeTets2, 'rows');
sharedTets = nodeTets1(tetIds, :);
valence = size(sharedTets, 1);

tetIds = find(ismember(sort(tets, 2), sharedTets, 'rows'));
end

