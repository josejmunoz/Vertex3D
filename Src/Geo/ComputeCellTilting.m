function [tiltingFaces] = ComputeCellTilting(Cell)
%COMPUTECELLTILTING Summary of this function goes here
%   Detailed explanation goes here
    tiltingFaces = [];
    for face = Cell.Faces
        for tris = face.Tris
            if length(tris.SharedByCells) > 2 && tris.Location == 1
                tiltingFaces = [tiltingFaces, ComputeEdgeTilting(tris, Cell.Y)];
            end
        end
    end
    if isempty(tiltingFaces)
        tiltingFaces = -1;
    end
end