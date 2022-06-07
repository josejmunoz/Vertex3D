function [features] = ComputeCellTriFeatures(cell, Set)
%COMPUTECELLTRIFEATURES Summary of this function goes here
%   Detailed explanation goes here
    features = struct();
    
    % Compute different measurements from the CELLS.Tris
    totalTris = 1;
    for face = cell.Faces
        [energyTris] = ComputeTriEnergy(face, cell.Y, Set);
        [~, areaTris] = ComputeFaceArea(vertcat(face.Tris.Edge), cell.Y, face.Centre);
        [~, perimeterTris] = ComputeFacePerimeter(vertcat(face.Tris.Edge), cell.Y, face.Centre);
        
        for numTris = 1:length(face.Tris)
            features(totalTris).energyTris = energyTris(numTris);
            features(totalTris).areaTris = areaTris{numTris};
            features(totalTris).perimeterTris = perimeterTris{numTris};
            features(totalTris).aspectRatioTris = areaTris{numTris} / perimeterTris{numTris};
            
            totalTris = totalTris + 1;
        end
        
    end
end

