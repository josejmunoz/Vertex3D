function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip23(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)	
	for c = 1:Geo.nCells
		for f = 1:length(Geo.Cells(c).Faces)
		    Ys = Geo.Cells(c).Y;
		    Ts = Geo.Cells(c).T;
			Face = Geo.Cells(c).Faces(f);
			nrgs = ComputeTriEnergy(Face, Ys, Set);
			Geo_backup = Geo; Geo_n_backup = Geo_n;
			for t = 1:length(Face.Tris)
				if ismember(Geo.Cells(c).globalIds(Face.Tris(t).Edge(1)),newYgIds)
					nrgs(t) = 0;
				end
			end
			
			[~,idVertex]=max(nrgs);
			YsToChange = Face.Tris(idVertex).Edge;
			
			if max(nrgs)<Set.RemodelTol || length(unique([Face.Tris.Edge])) == 3 || ...
					CheckSkinnyTriangles(Ys(YsToChange(1),:),Ys(YsToChange(2),:),Face.Centre)
                continue
			end
			
			[Ynew, Tnew] = YFlip23(Ys, Ts, YsToChange, Geo);

			ghostNodes = ismember(Tnew,Geo.XgID);
			ghostNodes = all(ghostNodes,2);
			if any(ghostNodes)
				fprintf('=>> Flips 2-2 are not allowed for now\n');
				return
			end

			targetTets = Geo.Cells(c).T(YsToChange,:);
			Geo   = ReplaceYs(targetTets, Tnew, Ynew, Geo);
			Geo_n = ReplaceYs(targetTets, Tnew, Ynew, Geo_n);
			
			Geo   = Rebuild(Geo, Set); 
			Geo_n = Rebuild(Geo_n, Set);
			
	        Geo   = BuildGlobalIds(Geo); 
			Geo_n = BuildGlobalIds(Geo_n);
			
			Geo   = UpdateMeasures(Geo);
			Geo_n = UpdateMeasures(Geo_n);

            if ~CheckConvexity(Tnew, Geo_backup) && CheckTris(Geo)
				fprintf('=>> 23 Flip.\n');
				Dofs = GetDOFs(Geo, Set);
				[Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
				[Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
				if DidNotConverge
					Geo   = Geo_backup;
					Geo_n = Geo_n_backup;
					fprintf('=>> 23-Flip rejected: did not converge\n');
					continue
				end
				newYgIds = unique([newYgIds; Geo.AssemblegIds]);
				Geo   = UpdateMeasures(Geo);
				Geo_n = UpdateMeasures(Geo_n);
			else
            	Geo   = Geo_backup;
				Geo_n = Geo_n_backup;
    		    fprintf('=>> 23-Flip rejected: is not compatible\n');
				continue
            end
		end
	end
end

