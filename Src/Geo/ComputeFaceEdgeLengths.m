function [EdgeLength, LengthsToCentre, AspectRatio] = ComputeFaceEdgeLengths(Face, Ys)
%% Compute the length of the edges of a face
    for currentTri = 1:length(Face.Tris)
        [EdgeLength{currentTri}, LengthsToCentre{currentTri}, AspectRatio{currentTri}] = ComputeTriLengthMeasurements(Face.Tris, Ys, currentTri, Face.Centre);
    end
end