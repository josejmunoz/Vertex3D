function Geo = BuildXFromY(Geo_n, Geo)
    for c = 1:length(Geo.Cells)
		% TODO FIXME, seems not optimal.. 2 loops necessary ?
        if ~isempty(Geo.Cells(c).T)
            if ismember(c, Geo.XgID) % Updating Ghost node
                dY = zeros(size(Geo.Cells(c).T,1), 3);
                for tet = 1:size(Geo.Cells(c).T,1)
                    gTet = Geo.Cells(c).T(tet,:);
                    for cm = 1:Geo.nCells
                        Cell   = Geo.Cells(cm);
                        Cell_n = Geo_n.Cells(cm);
                        hit = sum(ismember(Cell.T,gTet),2)==4;
                        if any(hit)
                            dY(tet,:) = Cell.Y(hit,:) - Cell_n.Y(hit,:);
                        end
                    end
                end
                Geo.Cells(c).X = Geo.Cells(c).X + mean(dY);
            else % Updating a main node
                dY = Geo.Cells(c).Y - Geo_n.Cells(c).Y;
                Geo.Cells(c).X = Geo.Cells(c).X + mean(dY);
            end
        end
    end
end

