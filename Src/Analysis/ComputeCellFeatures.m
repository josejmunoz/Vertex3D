function [features] = ComputeCellFeatures(cell)
%COMPUTECELLFEATURES Summary of this function goes here
%   Detailed explanation goes here
    features = struct();
    % Compute different measurements from the CELLS
    features.Area = ComputeCellArea(cell);
    features.Vol = ComputeCellVolume(cell);
    features.Height = ComputeCellHeight(cell);
    features.Area_Top = ComputeCellArea(cell, 'Top');
    features.Area_Bottom = ComputeCellArea(cell, 'Bottom');
    features.Area_CellCell = ComputeCellArea(cell, 'CellCell');
    features.Neighbours = length(ComputeCellNeighbours(cell));
    features.Neighbours_Top = length(ComputeCellNeighbours(cell, 'Top'));
    features.Neighbours_Bottom = length(ComputeCellNeighbours(cell, 'Bottom'));
    features.Tilting = mean(ComputeCellTilting(cell));

    %TODO: Other cell measurements
    %ComputeCellCircularity

end

