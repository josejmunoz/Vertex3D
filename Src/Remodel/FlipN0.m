function [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = FlipN0(Geo, Geo_n, Geo_0, Dofs, newYgIds, nodeToRemove, nodeToKeep, Set)
%FLIPN0 Summary of this function goes here
%   Detailed explanation goes here
flipName = 'Flip N-0';
oldTets = Geo.Cells(nodeToRemove).T;
nodesToCombine = [nodeToKeep, nodeToRemove];

oldYs = cellfun(@(x) GetYFromTet(Geo, x), num2cell(oldTets, 2), 'UniformOutput', false);
oldYs = vertcat(oldYs{:});

[newGeo, Tnew, Ynew] = CombineTwoGhostNodes(Geo, Set, nodesToCombine, oldTets, oldYs);
Geo.Cells(nodesToCombine(1)).X = newGeo.Cells(nodesToCombine(1)).X;

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);
end
end

