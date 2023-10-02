function [perimeter] = ComputePerimeterByLocation(cell, location)
%COMPUTEPERIMETERBYLOCATION Summary of this function goes here
%   Detailed explanation goes here
    perimeter = 0;
    for face = cell.Faces
        if face.InterfaceType == location
            for t = 1:length(face.Tris)
                currentTri = face.Tris(t);
                if length(currentTri.SharedByCells) > 1
                    perimeter = perimeter + currentTri.EdgeLength;
                end
            end
        end
    end
end

