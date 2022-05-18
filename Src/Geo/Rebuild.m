function Geo = Rebuild(Geo, Set)
	% TODO FIXME, whole function needs to be rethought
	for cc = 1:Geo.nCells
        Cell = Geo.Cells(cc);
        Neigh_nodes = unique(Geo.Cells(cc).T);
        Neigh_nodes(Neigh_nodes==cc)=[];
		for j  = 1:length(Neigh_nodes)
	        cj    = Neigh_nodes(j);
            ij			= [cc, cj];
            face_ids	= sum(ismember(Cell.T,ij),2)==2;
			% TODO FIXME BAD PROGRAMMING...
			newFace = true;
			for jj = 1:length(Geo.Cells(cc).Faces)
				Face = Geo.Cells(cc).Faces(jj);
				if ismember(cj, Face.ij)
					newFace = false;
					break
				end
			end
			if newFace
				Geo.Cells(cc).Faces(j+1:length(Geo.Cells(cc).Faces)+1)=Geo.Cells(cc).Faces(j:length(Geo.Cells(cc).Faces));
				Geo.Cells(cc).Faces(j)=BuildFace(cc, cj, Geo.nCells, Geo.Cells(cc), Geo.XgID, Set);
				Geo.Cells(cc).Faces(j).Centre = sum(Geo.Cells(cc).Y(face_ids,:),1)/sum(face_ids);
			else
				% TODO FIXME, I think this is an unnecessary call most of the time...
                [Face.Tris] = BuildEdges(Geo.Cells(cc).T, face_ids, Geo.Cells(cc).Faces(j).Centre, Geo.Cells(cc).Faces(j).InterfaceType, Geo.Cells(cc).X, Geo.Cells(cc).Y, 1:Geo.nCells);
            	Geo.Cells(cc).Faces(j).InterfaceType	= BuildInterfaceType(ij, Geo.XgID, Geo.XgTop, Geo.XgBottom);
                Geo.Cells(cc).Faces(j).ij = ij;
				[Geo.Cells(cc).Faces(j).Area] = ComputeFaceArea(vertcat(Geo.Cells(cc).Faces(j).Tris.Edge), Geo.Cells(cc).Y, Geo.Cells(cc).Faces(j).Centre);	
			end
		end
		Geo.Cells(cc).Area  = ComputeCellArea(Geo.Cells(cc));
        Geo.Cells(cc).Vol   = ComputeCellVolume(Geo.Cells(cc));
	end
end