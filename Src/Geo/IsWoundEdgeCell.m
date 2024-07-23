function [booleanWoundEdgeCell, booleanWoundEdgeCell_Top, booleanWoundEdgeCell_Bottom, booleanDebrisCell] = IsWoundEdgeCell(cell, debrisCells)
%ISWOUNDEDGECELL Summary of this function goes here
%   Detailed explanation goes here
    if ismember(cell.ID, debrisCells)
        booleanDebrisCell = 1;
        booleanWoundEdgeCell = 0;
        booleanWoundEdgeCell_Top = 0;
        booleanWoundEdgeCell_Bottom = 0;
    else
        booleanDebrisCell = 0;
        booleanWoundEdgeCell = any(ismember(ComputeCellNeighbours(cell), debrisCells));
        booleanWoundEdgeCell_Top = any(ismember(ComputeCellNeighbours(cell, 0), debrisCells));
        booleanWoundEdgeCell_Bottom = any(ismember(ComputeCellNeighbours(cell, 2), debrisCells));
    end
end

