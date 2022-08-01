function [surfaceArea] = ComputeTetSArea(newTets, Xs)
%COMPUTETETSAREA Summary of this function goes here
%   Detailed explanation goes here
allTris = nchoosek(newTets, 3);
surfaceArea = 0;
for tris = allTris'
    [area] = ComputeFaceArea(tris(1:2)', Xs, Xs(tris(3), :));
    surfaceArea = surfaceArea + area;
end
end

