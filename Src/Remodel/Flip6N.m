function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip6N(segmentToChange, Ts, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP5N Summary of this function goes here
%   Detailed explanation goes here

Ys = GetYFromTet(Geo, Ts);
[Ynew, Tnew] = YFlip6N(Ys, Ts, segmentToChange, Geo);

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName);
end
end