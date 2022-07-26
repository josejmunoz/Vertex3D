function [sideLengths] = ComputeTriSideLengths(Face, trisToChange, numCell, Geo)
%COMPUTETRISIDELENGTHS Summary of this function goes here
%   1st position is the edge not connecting to the face centre
edgeLenghts = zeros(1, length(Face.Tris));
for numTris = 1:length(Face.Tris)
    edgeLenghts(numTris) = ComputeEdgeLength(Face.Tris(numTris).Edge, Geo.Cells(numCell).Y);
end

lengthsToCentre = pdist2(Geo.Cells(numCell).Y(Face.Tris(trisToChange).Edge, :), Face.Centre);

sideLengths = [edgeLenghts(trisToChange), lengthsToCentre'];
end

