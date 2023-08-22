function Geo_new = CheckYsAndFacesHaveNotChanged(Geo, newTets, Geo_new)
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    aliveCells = nonDeadCells([Geo.Cells(nonDeadCells).AliveStatus] == 1);
    debrisCells = nonDeadCells([Geo.Cells(nonDeadCells).AliveStatus] == 0);
    for cellId = [aliveCells, debrisCells]
        if sum(any(ismember(newTets, Geo.XgBottom), 2)) > sum(any(ismember(newTets, Geo.XgTop), 2))
            tetsToCheck = ~any(ismember(Geo.Cells(cellId).T, Geo.XgBottom), 2);
            tetsToCheck_new = ~any(ismember(Geo_new.Cells(cellId).T, Geo.XgBottom), 2);
            interfaceType = 'Bottom';
        else
            tetsToCheck = ~any(ismember(Geo.Cells(cellId).T, Geo.XgTop), 2);
            tetsToCheck_new = ~any(ismember(Geo_new.Cells(cellId).T, Geo.XgTop), 2);
            interfaceType = 'Top';
        end
        %% Check that vertices that were untouched are not changing.
        assert(isequal(Geo.Cells(cellId).Y(tetsToCheck & any(ismember(Geo.Cells(cellId).T, Geo.XgID), 2), :), Geo_new.Cells(cellId).Y(tetsToCheck_new & any(ismember(Geo_new.Cells(cellId).T, Geo.XgID), 2), :)))
        
        %% Check that faces that were not changed, are not changing.
        for face = Geo.Cells(cellId).Faces
            if face.InterfaceType ~= interfaceType || ~ismember(cellId, newTets(:))
                idWithGeo_new = ismember(vertcat(Geo_new.Cells(cellId).Faces.ij), face.ij, 'rows');
                assert(sum(idWithGeo_new)==1)
                if ~isequal(Geo_new.Cells(cellId).Faces(idWithGeo_new).Centre, face.Centre)
                    Geo_new.Cells(cellId).Faces(idWithGeo_new).Centre = face.Centre;
                end
                
                assert(isequal(Geo_new.Cells(cellId).Faces(idWithGeo_new).Centre, face.Centre))
            end
        end
    end
end