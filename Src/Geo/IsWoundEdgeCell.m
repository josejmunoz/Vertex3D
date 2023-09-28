function [booleanWoundEdgeCell, booleanWoundEdgeCell_Top, booleanWoundEdgeCell_Bottom, booleanDebrisCell] = IsWoundEdgeCell(cell, debrisCells)
%ISWOUNDEDGECELL Summary of this function goes here
%   Detailed explanation goes here
    if ismember(cell.ID, debrisCells)
        booleanDebrisCell = 1;
        booleanWoundEdgeCell = 0;
    else
        booleanDebrisCell = 0;

    end
end

