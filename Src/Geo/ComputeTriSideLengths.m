function [sideLengths] = ComputeTriSideLengths(Face, trisToChange, Ys)
%COMPUTETRISIDELENGTHS Summary of this function goes here
%   1st position is the edge not connecting to the face centre
edgeLenghts = zeros(1, length(Face.Tris));
for numTris = 1:length(Face.Tris)
    edgeLenghts(numTris) = ComputeEdgeLength(Face.Tris(numTris).Edge, Ys);
end

lengthsToCentre(1) = norm(Ys(Face.Tris(trisToChange).Edge(1), :) - Face.Centre);
lengthsToCentre(2) = norm(Ys(Face.Tris(trisToChange).Edge(2), :) - Face.Centre);

sideLengths = [edgeLenghts(trisToChange), lengthsToCentre];
end

