function [Geo, Dofs] = UpdateDOFsCompress(Geo,Set)
	maxY = Geo.Cells(1).Y(1,2);
	for c = 1:Geo.nCells
		hit = find(Geo.Cells(c).Y(:,2)>maxY);
		if ~isempty(hit)
			maxY = max(Geo.Cells(c).Y(hit,2));
		end
		for f = 1:length(Geo.Cells(c).Faces)
			Face = Geo.Cells(c).Faces(f);
			if Geo.Cells(c).Faces(f).Centre(2)>maxY
				maxY = Geo.Cells(c).Faces(f).Centre(2);
			end
		end
	end
	Set.VPrescribed = maxY-Set.dx/((Set.TStopBC-Set.TStartBC)/Set.dt);
	Dofs = GetDOFs(Geo, Set);
	[dimP, numP] = ind2sub([3, Geo.numY+Geo.numF+Geo.nCells],Dofs.FixP);
	for c = 1:Geo.nCells
		prescYi  = ismember(Geo.Cells(c).globalIds, numP);
		Geo.Cells(c).Y(prescYi,dimP) = Set.VPrescribed;

		% TODO FIXME, I think this is proof that face global ids
		% should be in the cell struct and not the face struct
		for gn = 1:length(numP)
			for f = 1:length(Geo.Cells(c).Faces)
				Face = Geo.Cells(c).Faces(f);
				if numP(gn)==Face.globalIds
					Geo.Cells(c).Faces(f).Centre(dimP(gn)) = Set.VPrescribed;
				end
			end
		end
	end
end