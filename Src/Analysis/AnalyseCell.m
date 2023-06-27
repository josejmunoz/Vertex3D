function [c_features] = AnalyseCell(Geo, currentCell)
%ANALYSECELL Summary of this function goes here
%   Detailed explanation goes here
c_features = ComputeCellFeatures(Geo.Cells(currentCell));
c_features.ID = Geo.Cells(currentCell).ID;
if ismember(Geo.Cells(currentCell).ID, Geo.BorderCells)
    c_features.BorderCell = 1;
else
    c_features.BorderCell = 0;
end

end

