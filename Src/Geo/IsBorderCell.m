function boolBorderCell = IsBorderCell(Geo, currentCell)
boolBorderCell = ismember([Geo.Cells(currentCell).ID], Geo.BorderCells);
end