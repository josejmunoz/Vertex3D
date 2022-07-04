function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip32(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
	for c = 1:Geo.nCells
        f = 0;
		while f < length(Geo.Cells(c).Faces)
            f = f + 1;
    	    Ys = Geo.Cells(c).Y;
    	    Ts = Geo.Cells(c).T;

			Face = Geo.Cells(c).Faces(f);
			nrgs = ComputeTriEnergy(Face, Ys, Set);
			Geo_backup = Geo; Geo_n_backup = Geo_n;

            if max(nrgs)<Set.RemodelTol || length(unique([Face.Tris.Edge])) ~= 3 || ismember(Face.globalIds, newYgIds)
            	continue
            end

			YsToChange=[Face.Tris(1).Edge(1); Face.Tris(2).Edge(1); Face.Tris(3).Edge(1)];
            [Ynew, Tnew] = YFlip32(Ys, Ts, YsToChange, Geo);

            targetTets = Geo.Cells(c).T(YsToChange,:);
			Geo   = ReplaceYs(targetTets, Tnew, Ynew, Geo);
			Geo_n = ReplaceYs(targetTets, Tnew, Ynew, Geo_n);

            Geo   = RemoveFaces(f, Face.ij, Geo);
            Geo_n = RemoveFaces(f, Face.ij, Geo_n);

            %% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
			Geo   = Rebuild(Geo, Set);
			Geo_n = Rebuild(Geo_n, Set);

        	Geo   = BuildGlobalIds(Geo);
			Geo_n = BuildGlobalIds(Geo_n);

			Geo   = UpdateMeasures(Geo);
			Geo_n = UpdateMeasures(Geo_n);
            %% ----------------------------

            if ~CheckConvexity(Tnew,Geo_backup) && CheckTris(Geo)
    			fprintf('=>> 32 Flip.\n');
				Dofs = GetDOFs(Geo, Set);
				[Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
				[Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
				targetNodes = unique(targetTets);
				if DidNotConverge
					Geo   = Geo_backup;
					Geo_n = Geo_n_backup;
					fprintf('=>> 32-Flip rejected: did not converge\n');
					continue
				end
				for n_i = 1:length(unique(targetTets))
			        tNode = targetNodes(n_i);
			        news = find(sum(ismember(Tnew,tNode)==1,2));
			        if ~ismember(tNode, Geo.XgID)
				        Geo_n.Cells(tNode).Y(end-length(news)+1:end,:) = Geo.Cells(tNode).Y(end-length(news)+1:end,:);
			        end
				end
				newYgIds = unique([newYgIds; Geo.AssemblegIds]);
				Geo   = UpdateMeasures(Geo);
				Geo_n = UpdateMeasures(Geo_n);
%         	    return
            else
                Geo   = Geo_backup;
				Geo_n = Geo_n_backup;
                fprintf('=>> 32-Flip rejected: is not compatible\n');
    			continue
            end
		end
	end
end

