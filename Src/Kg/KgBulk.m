function [g,K,EnergyBulk]=KgBulk(Geo_0, Geo, Set) 
 
	[g, K] = initializeKg(Geo, Set); 
	 
	EnergyBulk=0; 
	errorInverted = []; 
 
	for c=1:Geo.nCells 
        if Geo.Remodelling
            if ~ismember(c,Geo.AssembleNodes)
                continue
            end
        end
        
        if isempty(Geo.Cells(c).AliveStatus) || Geo.Cells(c).AliveStatus ~= 1
            continue
        end
        
		ge=zeros(size(g, 1), 1); 
		cellNuclei  = Geo.Cells(c).X; 
		cellNuclei0 = Geo_0.Cells(c).X; 
		Ys  = Geo.Cells(c).Y; 
		Ys_0 = Geo_0.Cells(c).Y; 
		for f=1:length(Geo.Cells(c).Faces) 
			Tris = Geo.Cells(c).Faces(f).Tris; 
			for t=1:length(Tris) 
				y1   = Ys(Tris(t).Edge(1),:); 
				y1_0 = Ys_0(Tris(t).Edge(1),:); 
				y2   = Ys(Tris(t).Edge(2),:); 
				y2_0 = Ys_0(Tris(t).Edge(2),:);
                y3 = Geo.Cells(c).Faces(f).Centre;
                y3_0 = Geo_0.Cells(c).Faces(f).Centre;
                n3 = Geo.Cells(c).Faces(f).globalIds;
				currentTet     = [y1; y2; y3; cellNuclei]; 
				currentTet0    = [y1_0; y2_0; y3_0; cellNuclei0]; 
				currentTet_ids = [Geo.Cells(c).globalIds(Tris(t).Edge)', n3, Geo.Cells(c).cglobalIds]; 
				if Geo.Remodelling 
					if ~any(ismember(currentTet_ids,Geo.AssemblegIds)) 
                        continue
					end 
				end 
                try
                    [gB, KB, Energye] = KgBulkElem(currentTet, currentTet0, Set.mu_bulk, Set.lambda_bulk);
                    EnergyBulk=EnergyBulk+Energye;
                    ge=Assembleg(ge,gB,currentTet_ids);
                    K = AssembleK(K,KB,currentTet_ids);
                catch ME
                    if (strcmp(ME.identifier,'KgBulkElem:invertedTetrahedralElement'))
                        errorInverted = [errorInverted; currentTet_ids];
                    else
                        ME.rethrow();
                    end
                end
			end 
		end 
		g=g+ge; 
	end 
	if isempty(errorInverted) == 0 
		warning('Inverted Tetrahedral Element [%s]', sprintf('%d;', errorInverted')); 
	end 
end