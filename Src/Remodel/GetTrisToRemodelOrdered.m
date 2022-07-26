function [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here
energyPerCellAndFaces = [];
for c = 1:Geo.nCells
    Ys = Geo.Cells(c).Y;
    for numFace = 1:length(Geo.Cells(c).Faces)
        face = Geo.Cells(c).Faces(numFace);
        [nrgs]=ComputeTriEnergy(face, Ys, Set);
        if max(nrgs) >= Set.RemodelTol
            energyPerCellAndFaces(end+1, 1:3) = horzcat(c, numFace, max(nrgs));
        end
    end
end

if ~isempty(energyPerCellAndFaces)
    [energyPerCellAndFaces] = sortrows(energyPerCellAndFaces, 3, 'descend');
end

end

