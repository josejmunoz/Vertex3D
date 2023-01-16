function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip4N(segmentToChange, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP4N Summary of this function goes here
%   Detailed explanation goes here
hasConverged = false;
flipName = '4-4';
[Ynew, Tnew] = YFlip44(oldTets, segmentToChange, Geo, Set);

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName, segmentToChange);
end
end