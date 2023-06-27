function [nodeNeighbours] = getNodeNeighboursPerDomain(Geo, node, nodeOfDomain, mainNode)
%GETNODENEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here

allNodeTets = vertcat(Geo.Cells(node).T);

if ismember(nodeOfDomain, Geo.XgBottom)
    XgDomain = Geo.XgBottom;
elseif ismember(nodeOfDomain, Geo.XgTop)
    XgDomain = Geo.XgTop;
else
    XgDomain = Geo.XgLateral;
end

allNodeTets = allNodeTets(any(ismember(allNodeTets, XgDomain), 2), :);

if exist('mainNode', 'var')
    nodeNeighbours = unique(allNodeTets(any(ismember(allNodeTets, mainNode), 2), :));
else
    nodeNeighbours = unique(allNodeTets);
end

nodeNeighbours(ismember(nodeNeighbours, node)) = [];

end