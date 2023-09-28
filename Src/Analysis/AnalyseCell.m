function [c_features] = AnalyseCell(Geo, currentCell)
%ANALYSECELL Summary of this function goes here
%   Detailed explanation goes here
c_features = ComputeCellFeatures(Geo.Cells(currentCell));
c_features.ID = Geo.Cells(currentCell).ID;
c_features.BorderCell = IsBorderCell(Geo, currentCell);

% Compute different measurements from the WOUND
allCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus}))];
debrisCells = [allCells([allCells.AliveStatus] == 0).ID];

[c_features.WoundEdgeCell, c_features.WoundEdgeCell_Top, ...
    c_features.WoundEdgeCell_Bottom, c_features.DebrisCell] = IsWoundEdgeCell(cell, debrisCells);
end

