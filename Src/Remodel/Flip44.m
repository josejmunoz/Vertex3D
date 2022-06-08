function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip44(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
   for c = 1:Geo.nCells
		for f = 1:length(Geo.Cells(c).Faces)
	        Ys = Geo.Cells(c).Y;
	        Ts = Geo.Cells(c).T;
			% TODO FIXME, this should be a while? What happens with flip23?
			if f > length(Geo.Cells(c).Faces)
				break
			end
			Face = Geo.Cells(c).Faces(f);
			nrgs = ComputeTriEnergy(Face, Ys, Set);
			Geo_backup = Geo; Geo_n_backup = Geo_n;
		
			% The last ismember condition is necessary if a previous 44 flip
			% happens, as the new remodelled and approximated face might 
			% be in a cell not yet reviewed.
            % TODO: WHY THE RANGE OF REMODELING IS SET?
            % THAT MEANS THE 4-4 IS IN THE SAME FACE? SO NO 4-4 ARE ALLOWED
            % BETWEEN FACES?
			if max(nrgs)<Set.RemodelTol || min(nrgs)<Set.RemodelTol*1e-4 || length(unique([Face.Tris.Edge]))~=4 || ismember(Face.globalIds, newYgIds)
                continue
			end
			
			YsToChange=[Face.Tris(1).Edge(1); Face.Tris(2).Edge(1); Face.Tris(3).Edge(1); Face.Tris(4).Edge(1)];
            [Ynew, Tnew] = YFlip44(Ys, Ts, YsToChange, Face, Geo);
				
            targetTets = Geo.Cells(c).T(YsToChange,:);
			Geo   = ReplaceYs(targetTets, Tnew, Ynew, Geo);
			Geo_n = ReplaceYs(targetTets, Tnew, Ynew, Geo_n);

            Geo   = RemoveFaces(f, Face.ij, Geo);
            Geo_n = RemoveFaces(f, Face.ij, Geo_n);
			
			Geo   = Rebuild(Geo, Set); 
			Geo_n = Rebuild(Geo_n, Set);
			
	        Geo   = BuildGlobalIds(Geo); 
			Geo_n = BuildGlobalIds(Geo_n);

			Geo   = UpdateMeasures(Geo);
			Geo_n = UpdateMeasures(Geo_n);

            if ~CheckConvexity(Tnew,Geo_backup) && CheckTris(Geo)
    			fprintf('=>> 44 Flip.\n');
				Dofs = GetDOFs(Geo, Set);
				[Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
				[Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
				if DidNotConverge
					Geo   = Geo_backup;
					Geo_n = Geo_n_backup;
					fprintf('=>> 44-Flip rejected: did not converge\n');
					continue
				end
				newYgIds = unique([newYgIds; Geo.AssemblegIds]);
				targetNodes = unique(targetTets);
                for n_i = 1:length(targetNodes)
			        tNode = targetNodes(n_i);
			        news = find(sum(ismember(Tnew,tNode)==1,2));
			        if ~ismember(tNode, Geo.XgID)
				        Geo_n.Cells(tNode).Y(end-length(news)+1:end,:) = Geo.Cells(tNode).Y(end-length(news)+1:end,:);
			        end
                end
				Geo   = UpdateMeasures(Geo);
				Geo_n = UpdateMeasures(Geo_n);
            else
                Geo   = Geo_backup;
				Geo_n = Geo_n_backup;
                fprintf('=>> 44-Flip rejected: is not compatible\n');
    			continue
            end
		end
    end
end


