function [energyPerCellAndFaces] = GetTrisToRemodelOrdered(Geo, Set)
%GETTRISTOREMODELORDERED Summary of this function goes here
%   Detailed explanation goes here
energyPerCellAndFaces = [];
for c = 1:Geo.nCells
    for numFace = 1:length(Geo.Cells(c).Faces)
        face = Geo.Cells(c).Faces(numFace);
        %[nrgs]=ComputeTriEnergy(face, Ys, Set);
        for numTri = 1:length(face.Tris)
            [sideLengths] = ComputeTriSideLengths(face, numTri, Geo.Cells(c).Y);
            [aspectRatio(numTri)] = ComputeTriAspectRatio(sideLengths);
        end
        if max(aspectRatio) >= Set.RemodelTol
            energyPerCellAndFaces(end+1, 1:4) = horzcat(c, numFace, max(aspectRatio), face.globalIds);
        end
    end
end

if ~isempty(energyPerCellAndFaces)
    [energyPerCellAndFaces] = sortrows(energyPerCellAndFaces, 3, 'descend');
end

end

