function Geo = BuildXFromY(Geo_n, Geo)

    proportionOfMax = 0;

    aliveCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    allCellsToUpdate = setdiff(1:length(Geo.Cells), vertcat(Geo.BorderCells, Geo.BorderGhostNodes));
    for c = allCellsToUpdate
		% TODO FIXME, seems not optimal.. 2 loops necessary ?
        if ~isempty(Geo.Cells(c).T)
            if ismember(c, Geo.XgID) % Updating Ghost node
                dY = zeros(size(Geo.Cells(c).T,1), 3);
                for tet = 1:size(Geo.Cells(c).T,1)
                    gTet = Geo.Cells(c).T(tet,:);
                    gTet_Cells = gTet(ismember(gTet, aliveCells));
                    cm = gTet_Cells(1);
                    Cell   = Geo.Cells(cm);
                    Cell_n = Geo_n.Cells(cm);
                    hit = sum(ismember(Cell.T,gTet),2)==4;
                    dY(tet,:) = Cell.Y(hit,:) - Cell_n.Y(hit,:);
                end
                Geo.Cells(c).X = Geo.Cells(c).X + (proportionOfMax) * max(dY) + (1 - proportionOfMax) * mean(dY);
            else % Updating a main node
                dY = Geo.Cells(c).Y - Geo_n.Cells(c).Y;
                Geo.Cells(c).X = Geo.Cells(c).X + (proportionOfMax) * max(dY) + (1 - proportionOfMax) * mean(dY);
            end
        end
    end
end

