function [newTets] = DoFlip13(Geo, newNodeIDs, oldTets)
%DOFLIP13 Summary of this function goes here
%   Detailed explanation goes here

mainNode = oldTets(~cellfun(@isempty, {Geo.Cells(oldTets).AliveStatus}));
remainingNodes = setdiff(oldTets, mainNode);

if length(mainNode) > 1
    error('DoFlip13 more than 1 mainNode');
end

combinationsToCreate = nchoosek(remainingNodes, 2);

newTets = [];
for combination = combinationsToCreate'
    newTets(end+1, 1:4) = [newNodeIDs, mainNode, combination'];
end
end

