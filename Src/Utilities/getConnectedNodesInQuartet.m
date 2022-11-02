function [connectedNodes, unconnectedNodes] = getConnectedNodesInQuartet(Geo, Xs)
%GETCONNECTEDNODES Summary of this function goes here
%   Detailed explanation goes here
cellNode_Adjacency = cell(length(Xs), 1);
for numNode = 1:length(Xs)
    currentNeighbours = getNodeNeighbours(Geo, Xs(numNode));
    cellNode_Adjacency{numNode} = Xs(ismember(Xs, [currentNeighbours; Xs(numNode)]));
end

connectedNodes = Xs(cellfun(@length, cellNode_Adjacency) == length(Xs));
unconnectedNodes = Xs(cellfun(@length, cellNode_Adjacency) < length(Xs));
end

