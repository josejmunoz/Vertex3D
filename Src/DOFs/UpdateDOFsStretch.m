function Geo = UpdateDOFsStretch(FixP, Geo, Set)
	for c = 1:Geo.nCells
		prescYi  = ismember(Geo.Cells(c).globalIds, FixP);
		Geo.Cells(c).Y(prescYi,2) = Geo.Cells(c).Y(prescYi,2) + Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
		% TODO FIXME, I think this is proof that face global ids
		% should be in the cell struct and not the face struct
		for gn = 1:length(FixP)
			for f = 1:length(Geo.Cells(c).Faces)
				Face = Geo.Cells(c).Faces(f);
				if length(Face.Tris)==3
					continue
				end
				if FixP(gn)==Face.globalIds
					Geo.Cells(c).Faces(f).Centre(2) = Geo.Cells(c).Faces(f).Centre(2) + Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
				end
			end
		end
	end
end