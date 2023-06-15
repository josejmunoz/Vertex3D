function Geo = RemoveFaces(f, ij, Geo)
    oppfaceId = [];
    for f2 = 1:length(Geo.Cells(ij(2)).Faces)
	    Faces2 = Geo.Cells(ij(2)).Faces(f2);
	    if all(ismember(ij, Faces2.ij))
		    oppfaceId = f2;
	    end
    end
    Geo.Cells(ij(1)).Faces(f) = [];
    if ~isempty(oppfaceId)
	    Geo.Cells(ij(2)).Faces(oppfaceId) = [];
    end
end