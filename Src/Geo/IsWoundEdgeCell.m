function [booleanWoundEdgeCell, booleanWoundEdgeCell_Top, booleanWoundEdgeCell_Bottom, booleanDebrisCell] = IsWoundEdgeCell(cell, debrisCells)
%ISWOUNDEDGECELL Summary of this function goes here
%   Detailed explanation goes here
    if ismember(cell.ID, debrisCells)
        booleanDebrisCell = 1;
        booleanWoundEdgeCell = 0;
    else
        booleanDebrisCell = 0;
        booleanWoundEdgeCell = any(ismember(ComputeCellNeighbours(cell), debrisCells));
        booleanWoundEdgeCell_Top = any(ismember(ComputeCellNeighbours(cell, "Top"), debrisCells));
        booleanWoundEdgeCell_Bottom = any(ismember(ComputeCellNeighbours(cell, "Bottom"), debrisCells));
    end
end

