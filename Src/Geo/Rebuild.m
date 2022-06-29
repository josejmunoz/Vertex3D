function Geo = Rebuild(Geo, Set)
	% TODO FIXME, whole function needs to be rethought
    oldGeo = Geo;
	for cc = 1:Geo.nCells
        Cell = Geo.Cells(cc);
        Neigh_nodes = unique(Geo.Cells(cc).T);
        Neigh_nodes(Neigh_nodes==cc)=[];
		for j  = 1:length(Neigh_nodes)
	        cj    = Neigh_nodes(j);
            ij			= [cc, cj];
            face_ids	= sum(ismember(Cell.T,ij),2)==2;
            [oldFaceExists, previousFace] = ismember(cj, [oldGeo.Cells(cc).Faces.ij]);
            
			Geo.Cells(cc).Faces(j) = BuildFace(cc, cj, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set, Geo.XgTop, Geo.XgBottom);
            if ~oldFaceExists
                Geo.Cells(cc).Faces(j).Centre = BuildFaceCentre(ij, Geo.nCells, Geo.Cells(cc).X, Geo.Cells(cc).Y(face_ids,:), Set.f);
            else
                previousFace = ceil(previousFace/2);
                Geo.Cells(cc).Faces(j).Centre = oldGeo.Cells(cc).Faces(previousFace).Centre;
            end
        end
        Geo.Cells(cc).Faces = Geo.Cells(cc).Faces(1:length(Neigh_nodes));
		Geo.Cells(cc).Area  = ComputeCellArea(Geo.Cells(cc));
        Geo.Cells(cc).Vol   = ComputeCellVolume(Geo.Cells(cc));
	end
end