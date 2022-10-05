function [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here
energyPerCellAndFaces = [];
for c = 1:Geo.nCells
    for numFace = 1:length(Geo.Cells(c).Faces)
        face = Geo.Cells(c).Faces(numFace);
        
        if median([face.Tris.AspectRatio]) >= Set.RemodelTol && ~isequal(face.InterfaceType, 'CellCell') 
            energyPerCellAndFaces(end+1, 1:6) = horzcat(c, numFace, max([face.Tris.AspectRatio]), face.globalIds, face.ij);
        end
    end
end

if ~isempty(energyPerCellAndFaces)
    [energyPerCellAndFaces] = sortrows(energyPerCellAndFaces, 3, 'descend');
end

end

