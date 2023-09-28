function boolBorderCell = IsBorderCell(Geo, currentCell)
if ismember(Geo.Cells(currentCell).ID, Geo.BorderCells)
    boolBorderCell = 1;
else
    boolBorderCell = 0;
end
end