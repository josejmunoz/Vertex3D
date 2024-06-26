function [g,K,Energy_T]=KgSurfaceCellBasedAdhesion(Geo, Set)
	[g, K] = initializeKg(Geo, Set);
	Energy_T = 0;
	for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
		if Geo.Remodelling
			if ~ismember(c,Geo.AssembleNodes)
        		continue
			end
		end
%         if Geo.Cells(c).AliveStatus ~= 1
%             continue
%         end
        
        Energy_c = 0;

		Cell  = Geo.Cells(c);
		Ys    = Geo.Cells(c).Y;
		ge	  = sparse(size(g, 1), 1);
		fact0 = 0;
		for f=1:length(Cell.Faces)
			face = Cell.Faces(f);
			if face.InterfaceType == 1
				Lambda=Set.lambdaS1*Cell.ExternalLambda;
			elseif face.InterfaceType == 2
				Lambda=Set.lambdaS2*Cell.InternalLambda;
			elseif face.InterfaceType == 3
				Lambda=Set.lambdaS3*Cell.SubstrateLambda;
			end
			fact0=fact0+Lambda*face.Area;
		end
		fact=fact0/Cell.Area0^2;
        for f=1:length(Cell.Faces)
			face = Cell.Faces(f);
			Tris=Cell.Faces(f).Tris;
			if face.InterfaceType == 1
				Lambda=Set.lambdaS1*Cell.ExternalLambda;
			elseif face.InterfaceType == 2
				Lambda=Set.lambdaS2*Cell.InternalLambda;
			elseif face.InterfaceType == 3
				Lambda=Set.lambdaS3*Cell.SubstrateLambda;
			end
            for t = 1:length(Tris)
				y1 = Ys(Tris(t).Edge(1),:);
				y2 = Ys(Tris(t).Edge(2),:);
				y3 = Cell.Faces(f).Centre;
				n3 = Cell.Faces(f).globalIds;
				nY = [Cell.globalIds(Tris(t).Edge)', n3];
				if Geo.Remodelling
					if ~any(ismember(nY,Geo.AssemblegIds))
                        continue
					end
				end
				[gs,Ks,Kss]=gKSArea(y1,y2,y3);
				gs=Lambda*gs;
            	ge=Assembleg(ge,gs,nY);
				Ks=fact*Lambda*(Ks+Kss);
				K = AssembleK(K,Ks,nY);
            end
        end
		g=g+ge*fact;
		K=K+(ge)*(ge')/(Cell.Area0^2);
    	Energy_c=Energy_c+ (1/2)*fact0*fact;
        Energy(c) = Energy_c;
    end
    Energy_T = sum(Energy);
end

