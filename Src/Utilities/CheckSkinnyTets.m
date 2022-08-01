function [skinnyTets, aspectRatio] = CheckSkinnyTets(newTets, Geo)
%CHECKSKINNYTETS Summary of this function goes here
%   Detailed explanation goes here

aspectRatio = [];
for tet = newTets'
    [SArea] = ComputeTetSArea(tet, vertcat(Geo.Cells.X));
    [vol] = ComputeTetVolume(tet, Geo);
    aspectRatio(end+1) = SArea/vol;
end

skinnyTets = aspectRatio > 70;

end

