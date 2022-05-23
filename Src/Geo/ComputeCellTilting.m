function [tiltingFaces] = ComputeCellTilting(Cell)
%COMPUTECELLTILTING Summary of this function goes here
%   Detailed explanation goes here
    tiltingFaces = [];
    for face = Cell.Faces
        for tris = face.Tris
            if length(tris.SharedByCells) > 2 && tris.Location == 'Cell-Cell'
                v1 = Cell.Y(tris.Edge(1), :) - Cell.Y(tris.Edge(2), :); % realEdge
                fixedVertex = [Cell.Y(tris.Edge(1), 1:2), Cell.Y(tris.Edge(2), 3)];
                %TODO: CHECK IF THIS IS CORRECT
                %TODO: Improve perpendicular edge for curve tissues
                v2 = Cell.Y(tris.Edge(1), :) - fixedVertex;% Perpendicular edge
                % Calculate the angle between the perpendicular edge and
                % the real one.
                tiltingFaces = [tiltingFaces, atan2(norm(cross(v1,v2)),dot(v1,v2)) * 100];
            end
        end
    end
end