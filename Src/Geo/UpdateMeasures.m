function Geo = UpdateMeasures(Geo)
	for c = 1:Geo.nCells
		for f = 1:length(Geo.Cells(c).Faces)
            [Geo.Cells(c).Faces(f).Area, triAreas]  = ComputeFaceArea(vertcat(Geo.Cells(c).Faces(f).Tris.Edge), Geo.Cells(c).Y, Geo.Cells(c).Faces(f).Centre);
            [Geo.Cells(c).Faces(f).Tris.Area] = triAreas{:};
            [edgeLengths] = ComputeFaceEdgeLengths(Geo.Cells(c).Faces(f), Geo.Cells(c).Y);
            [Geo.Cells(c).Faces(f).Tris.EdgeLength] = edgeLengths{:};
		end
		Geo.Cells(c).Area  = ComputeCellArea(Geo.Cells(c));
    	Geo.Cells(c).Vol   = ComputeCellVolume(Geo.Cells(c));
end