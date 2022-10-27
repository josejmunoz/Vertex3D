function [energiesPerCellAndFaces, allEnergies] = ComputeCellTriEnergy(Geo, Set)
%COMPUTECELLTRIENERGY Summary of this function goes here
%   Detailed explanation goes here

energiesPerCellAndFaces = table();
allEnergies = {};
for c = 1:Geo.nCells
    Ys = Geo.Cells(c).Y;
    for numFace = 1:length(Geo.Cells(c).Faces)
        face = Geo.Cells(c).Faces(numFace);
        [nrgs]=ComputeTriEnergy(face, Ys, Set);
        energiesPerCellAndFaces = vertcat(energiesPerCellAndFaces, table(c, numFace, max(nrgs)));
        allEnergies(end+1) = {nrgs};
    end
end

