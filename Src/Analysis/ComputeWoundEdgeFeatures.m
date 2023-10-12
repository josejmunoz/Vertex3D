function [woundEdgeFeatures] = ComputeWoundEdgeFeatures(Geo, debrisCells)
%COMPUTEWOUNDEDGEFEATURES Summary of this function goes here
%   Detailed explanation goes here
    %% Compute features
    if ~exist('debrisCells', 'var')
        debrisCells = getDebrisCells(Geo);
    end
    booleanWoundEdgeCell = [];
    for cell = Geo.Cells
        [booleanWoundEdgeCell(end+1)] = IsWoundEdgeCell(cell, debrisCells);
    end

    woundEdgeCells = Geo.Cells(booleanWoundEdgeCell == 1);
    woundEdgeFeatures = {};
    for woundEdgeCell = woundEdgeCells
        woundEdgeFeatures{end+1} = ComputeCellFeatures(woundEdgeCell);
    end
    woundEdgeFeatures = struct2table(vertcat(woundEdgeFeatures{:}));
    woundEdgeFeatures_mean = mean(table2array(woundEdgeFeatures));
    woundEdgeFeatures_mean = array2table(woundEdgeFeatures_mean, 'VariableNames', woundEdgeFeatures.Properties.VariableNames);
end

