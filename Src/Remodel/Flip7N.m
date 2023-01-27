function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged, Tnew] = Flip7N(segmentToChange, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP5N Summary of this function goes here
%   Detailed explanation goes here
hasConverged = false;
flipName = '7-N';
[Ynew, Tnew] = YFlip7N(oldTets, segmentToChange, Geo, Set);

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);
end
end