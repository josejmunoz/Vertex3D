function Geo = BuildXFromY(Geo_n, Geo)

    for c = 1:length(Geo.Cells)
		% TODO FIXME, seems not optimal.. 2 loops necessary ?
        if ~isempty(Geo.Cells(c).AliveStatus)
            if Geo.Cells(c).AliveStatus == 0 % Ghost node
                % Updating a ghost node
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
            else % Regular node
                % Updating a main node
                dY = Geo.Cells(c).Y - Geo_n.Cells(c).Y;
                Geo.Cells(c).X = Geo.Cells(c).X + mean(dY);
            end
        end
    end



%     for c = 1:length(Geo.Cells)
%         Ts = Geo.Cells(c).T;
%         currentTetrahedra = any(ismember(Ts, c), 2);
%         if any(currentTetrahedra)
%             if any(Geo.XgID==c)
%                 dYs = zeros(size(Geo.Cells(c).T,1),3);
%                 count = 1;
%                 for cm = 1:Geo.nCells
%                     gY_Id = sum(ismember(Geo.Cells(c).T, cm),2)==1;
%                     if any(gY_Id)
%                         dYs(count,:) = Geo.Cells(cm).Y(gY_Id,:) - Geo_n.Cells(cm).Y(gY_Id,:);
%                     end
%                 end
%             else
%                 dY = Geo.Cells(c).Y - Geo_n.Cells(c).Y;
%             end
%             changeOfSurroundingYs = mean(dY(currentTetrahedra, :));
%             Geo.Cells(c).X = Geo.Cells(c).X + changeOfSurroundingYs;
%         end
%     end
end

