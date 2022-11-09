function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip5N(segmentToChange, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP5N Summary of this function goes here
%   Detailed explanation goes here

flipName = '5-5';
[Ynew, Tnew] = YFlip55(oldTets, segmentToChange, Geo, Set);

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName, segmentToChange);
end
end