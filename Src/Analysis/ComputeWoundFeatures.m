function [features] = ComputeWoundFeatures(Geo, debrisCells)
%COMPUTEWOUNDFEATURES Summary of this function goes here
%   Detailed explanation goes here
    %% Init features
    features = struct();
    
    %% Compute features
    % Wound area: top and bottom
    if ~exist('debrisCells', 'var')
        debrisCells = getDebrisCells(Geo);
    end
    booleanWoundEdgeCell = [];
    for cell = Geo.Cells
        [booleanWoundEdgeCell(end+1)] = IsWoundEdgeCell(cell, debrisCells);
    end
    
    woundEdgeCells = Geo.Cells(booleanWoundEdgeCell == 1);
    borderVertices_Top = [];
    borderVertices_Bottom = [];
    for woundEdgeCell = woundEdgeCells
        for face = woundEdgeCell.Faces
            for tri = face.Tris
                if any(ismember(tri.SharedByCells, debrisCells))
                    if face.InterfaceType == 'Top'
                        borderVertices_Top = vertcat(borderVertices_Top, vertcat(woundEdgeCell.Y(tri.Edge, :)));
                    elseif face.InterfaceType == 'Bottom'
                        borderVertices_Bottom = vertcat(borderVertices_Bottom, vertcat(woundEdgeCell.Y(tri.Edge, :)));
                    end
                end
            end
        end
    end

    [features] = ComputeWoundEdgeFeatures(Geo, debrisCells);

    k = boundary(borderVertices_Top(:, 1), borderVertices_Top(:, 2), 0);
    features.woundArea_Top = polyarea(borderVertices_Top(k(1:end-1), 1), borderVertices_Top(k(1:end-1), 2));
    k = boundary(borderVertices_Bottom(:, 1), borderVertices_Bottom(:, 2), 0);
    features.woundArea_Bottom = polyarea(borderVertices_Bottom(k(1:end-1), 1), borderVertices_Bottom(k(1:end-1), 2));
    
end

