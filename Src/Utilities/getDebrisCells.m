function debrisCells = getDebrisCells(Geo)
    allCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus}))];
    debrisCells = [allCells([allCells.AliveStatus] == 0).ID];
end