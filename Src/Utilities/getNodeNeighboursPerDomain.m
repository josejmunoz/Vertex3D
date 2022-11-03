function [nodeNeighbours] = getNodeNeighboursPerDomain(Geo, node, domain, mainNode)
%GETNODENEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here

allNodeTets = vertcat(Geo.Cells(node).T);

if domain
    XgDomain = Geo.XgBottom;
else
    XgDomain = Geo.XgTop;
end

allNodeTets = allNodeTets(any(ismember(allNodeTets, XgDomain), 2), :);


if exist('mainNode', 'var')
    nodeNeighbours = unique(allNodeTets(any(ismember(allNodeTets, mainNode), 2), :));
else
    nodeNeighbours = unique(allNodeTets);
end

nodeNeighbours(ismember(nodeNeighbours, node)) = [];

end