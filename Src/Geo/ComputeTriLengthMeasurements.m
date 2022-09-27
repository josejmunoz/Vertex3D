function [EdgeLength, LengthsToCentre, AspectRatio] = ComputeTriLengthMeasurements(Tris, Ys, currentTri, FaceCentre)
%COMPUTETRILENGTHMEASUREMENTS Summary of this function goes here
%   Detailed explanation goes here
    EdgeLength = norm(Ys(Tris(currentTri).Edge(1), :) - Ys(Tris(currentTri).Edge(2), :));
    LengthsToCentre = [norm(Ys(Tris(currentTri).Edge(1), :) - FaceCentre), norm(Ys(Tris(currentTri).Edge(2), :) - FaceCentre)];
    AspectRatio = ComputeTriAspectRatio([EdgeLength LengthsToCentre]);
end

