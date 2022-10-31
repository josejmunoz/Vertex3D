function [Geo, Tnew, Ynew] = FlipN0(Geo, nodeToRemove, nodeToKeep, oldTets, Set)
%FLIPN0 Summary of this function goes here
%   Detailed explanation goes here
nodesToCombine = [nodeToKeep, nodeToRemove];
oldYs = cellfun(@(x) GetYFromTet(Geo, x), num2cell(oldTets, 2), 'UniformOutput', false);
oldYs = vertcat(oldYs{:});
[newGeo, Tnew, Ynew] = CombineTwoGhostNodes(Geo, Set, nodesToCombine, oldTets, oldYs);
Geo.Cells(nodesToCombine(1)).X = newGeo.Cells(nodesToCombine(1)).X;
end

