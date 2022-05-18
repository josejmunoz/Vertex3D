function Geo = UpdateVertices(Geo, Set, dy_reshaped)
    for c = 1:Geo.nCells
        dY = dy_reshaped(Geo.Cells(c).globalIds,:);
        Geo.Cells(c).Y = Geo.Cells(c).Y + dY;
		dYc = dy_reshaped(Geo.Cells(c).cglobalIds,:); 
		Geo.Cells(c).X = Geo.Cells(c).X + dYc; 
        for f = 1:length(Geo.Cells(c).Faces)
            Geo.Cells(c).Faces(f).Centre = Geo.Cells(c).Faces(f).Centre + dy_reshaped(Geo.Cells(c).Faces(f).globalIds,:);
        end
    end
end

