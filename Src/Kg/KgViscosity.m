function [g,K,EnergyF]=KgViscosity(Geo_n, Geo, Set)
    K=(Set.nu/Set.dt).*eye((Geo.numY+Geo.numF+Geo.nCells)*3);
	% TODO FIXME placeholder...
	dY = zeros(Geo.numF+Geo.numY+Geo.nCells,3);
	% TODO FIXME BAD!
	for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        % THERE WAS A HARD TO DEBUG ERROR HERE... 
		if Geo.Remodelling
			if ~ismember(c,Geo.AssembleNodes)
        		continue
			end
		end
		Cell = Geo.Cells(c);
		Cell_n = Geo_n.Cells(c);
		dY(Cell.globalIds,:) = (Cell.Y-Cell_n.Y);
		for f = 1:length(Cell.Faces)
			Face = Cell.Faces(f);
			Face_n = Cell_n.Faces(f);
            if ~isstring(Face.Centre) && ~isstring(Face_n.Centre)
			    dY(Face.globalIds,:) = (Face.Centre-Face_n.Centre);
            end
		end
		dY(Cell.cglobalIds,:) = (Cell.X-Cell_n.X);
	end
	g = (Set.nu/Set.dt).*reshape(dY', (Geo.numF+Geo.numY+Geo.nCells)*3, 1);
	EnergyF = (1/2)*(g')*g/Set.nu;

end