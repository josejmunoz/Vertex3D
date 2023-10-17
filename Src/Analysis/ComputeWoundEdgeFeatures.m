function [woundEdgeFeatures_mean] = ComputeWoundEdgeFeatures(Geo, woundEdgeCells)
%COMPUTEWOUNDEDGEFEATURES Summary of this function goes here
%   Detailed explanation goes here
    %% Compute features
    if ~exist('woundEdgeCells', 'var')
        woundEdgeCells = getDebrisCells(Geo);
        booleanWoundEdgeCell = [];
        for cell = Geo.Cells
            [booleanWoundEdgeCell(end+1)] = IsWoundEdgeCell(cell, woundEdgeCells);
        end
    else
        booleanWoundEdgeCell = ismember([Geo.Cells.ID], woundEdgeCells);
    end

    woundEdgeCells = Geo.Cells(booleanWoundEdgeCell == 1);
    woundEdgeFeatures = {};
    for woundEdgeCell = woundEdgeCells
        woundEdgeFeatures{end+1} = ComputeCellFeatures(woundEdgeCell);
    end
    woundEdgeFeatures = struct2table(vertcat(woundEdgeFeatures{:}));
    woundEdgeFeatures_mean = mean(table2array(woundEdgeFeatures));
    woundEdgeFeatures_mean = table2struct(array2table(woundEdgeFeatures_mean, 'VariableNames', woundEdgeFeatures.Properties.VariableNames));
    woundEdgeFeatures_mean.numberOfCells = length(woundEdgeCells);
end

