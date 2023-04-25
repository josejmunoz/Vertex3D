function [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged, Tnew] = FlipNM(segmentToChange, cellToIntercalateWith, oldTets, oldYs, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPNM Summary of this function goes here
%   Detailed explanation goes here
hasConverged = false;
flipName = 'N-M';
[Ynew, Tnew] = YFlipNM(oldTets, cellToIntercalateWith, oldYs, segmentToChange, Geo);

if ~isempty(Tnew)
    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName, segmentToChange);
end
end