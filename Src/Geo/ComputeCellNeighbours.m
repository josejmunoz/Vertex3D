function [neighbours] = ComputeCellNeighbours(cell, locationFilter)
%COMPUTECELLNEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here
    allSharedByCells = [];
    for face = cell.Faces
        if exist('filtering', 'var')
            if face.InterfaceType == locationFilter
                allSharedByCells = [allSharedByCells, face.Tris.SharedByCells];
            end
        else
            allSharedByCells = [allSharedByCells, face.Tris.SharedByCells];
        end
    end
    neighbours = unique(allSharedByCells);
    neighbours(neighbours == cell.ID) = [];
end

