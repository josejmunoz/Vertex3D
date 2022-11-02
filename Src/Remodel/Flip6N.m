function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip6N(segmentToChange, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP5N Summary of this function goes here
%   Detailed explanation goes here

flipName = 'Flip 6-N';

Ys = GetYFromTet(Geo, oldTets);
[Ynew, Tnew] = YFlip6N(Ys, oldTets, segmentToChange, Geo, Set);

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);
end
end