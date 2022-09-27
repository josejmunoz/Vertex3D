function [lengthsToCentre] = ComputeTriLengthsToCentre(Face, trisToChange, Ys)
%COMPUTETRILENGTHSTOCENTRE Summary of this function goes here
%   1st position is the edge not connecting to the face centre

lengthsToCentre(1) = norm(Ys(Face.Tris(trisToChange).Edge(1), :) - Face.Centre);
lengthsToCentre(2) = norm(Ys(Face.Tris(trisToChange).Edge(2), :) - Face.Centre);
end

