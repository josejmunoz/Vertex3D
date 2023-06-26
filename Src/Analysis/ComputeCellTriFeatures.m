function [features] = ComputeCellTriFeatures(cell, Set)
%COMPUTECELLTRIFEATURES Summary of this function goes here
%   Detailed explanation goes here
    features = struct();
    
    % Compute different measurements from the CELLS.Tris
    totalTris = 1;
    for face = cell.Faces
        [energyAreaTris] = ComputeTriAreaEnergy(face, Set);
        [energyARTris] = ComputeTriAREnergy(face, cell.Y, Set);
        [~, areaTris] = ComputeFaceArea(vertcat(face.Tris.Edge), cell.Y, face.Centre);
        [~, perimeterTris] = ComputeFacePerimeter(vertcat(face.Tris.Edge), cell.Y, face.Centre);
        
        
        for numTris = 1:length(face.Tris)
            features(totalTris).energyAreaTris = energyAreaTris(numTris);
            features(totalTris).energyARTris = energyARTris(numTris);
            features(totalTris).areaTris = areaTris{numTris};
            features(totalTris).perimeterTris = perimeterTris{numTris};
            features(totalTris).aspectRatioTris = face.Tris(numTris).AspectRatio;
            
            features(totalTris).facingNode = setdiff(face.ij, cell.ID);
            features(totalTris).faceType = face.InterfaceType;
            
            totalTris = totalTris + 1;
        end
        
    end
end

