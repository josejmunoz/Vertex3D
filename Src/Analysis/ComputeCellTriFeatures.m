function [features] = ComputeCellTriFeatures(cell, Set)
%COMPUTECELLTRIFEATURES Summary of this function goes here
%   Detailed explanation goes here
    features = struct();
    
    % Compute different measurements from the CELLS.Tris
    numFace = 1;
    energyTris = {};
    trisArea = {};
    trisPerimeter = {};
    for face = cell.Faces
        [energyTris{numFace}] = ComputeTriEnergy(face, cell.Y, Set);
        [~, trisArea{numFace}] = ComputeFaceArea(vertcat(face.Tris.Edge), cell.Y, face.Centre);
        [~, trisPerimeter{numFace}] = ComputeFacePerimeter(vertcat(face.Tris.Edge), cell.Y, face.Centre);
        numFace = numFace + 1;
    end
    features.energyTris = [energyTris{:}];
    features.areaTris = [trisArea{:}];
    features.perimeterTris = [trisPerimeter{:}];
end

