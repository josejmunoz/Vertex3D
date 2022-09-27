function [Tris] = ComputeTriLengthMeasurements(Tris, Ys, currentTri, FaceCentre)
%COMPUTETRILENGTHMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here
    Tris(currentTri).EdgeLength = norm(Ys(Tris(currentTri).Edge(1), :) - Ys(Tris(currentTri).Edge(2), :));
    Tris(currentTri).LengthsToCentre = [norm(Ys(Tris(currentTri).Edge(1), :) - FaceCentre), norm(Ys(Tris(currentTri).Edge(2), :) - FaceCentre)];
    Tris(currentTri).AspectRatio = ComputeTriAspectRatio([Tris(currentTri).EdgeLength Tris(currentTri).LengthsToCentre]);
end

