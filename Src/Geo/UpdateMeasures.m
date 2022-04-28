function Geo = UpdateMeasures(Geo)
	for c = 1:Geo.nCells
		for f = 1:length(Geo.Cells(c).Faces)
        	[Geo.Cells(c).Faces(f).Area, Geo.Cells(c).Faces(f).TrisArea] = ComputeFaceArea(Geo.Cells(c).Faces(f), Geo.Cells(c).Y);
		end
		Geo.Cells(c).Area  = ComputeCellArea(Geo.Cells(c));
    	Geo.Cells(c).Vol   = ComputeCellVolume(Geo.Cells(c));
end