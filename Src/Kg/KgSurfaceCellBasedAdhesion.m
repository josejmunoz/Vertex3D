function [g,K,EnergyS]=KgSurfaceCellBasedAdhesion(Geo, Set)
	[g, K] = initializeKg(Geo, Set);
	EnergyS = 0;
	for c = 1:Geo.nCells
		if Geo.Remodelling
			if ~ismember(c,Geo.AssembleNodes)
        		continue
			end
		end
        if isempty(Geo.Cells(c).AliveStatus) || Geo.Cells(c).AliveStatus ~= 1
            continue
        end
        
		Cell  = Geo.Cells(c);
		Ys    = Geo.Cells(c).Y;
		ge	  = sparse(size(g, 1), 1);
		fact0 = 0;
		for f=1:length(Cell.Faces)
			face = Cell.Faces(f);
			if face.InterfaceType == 'Top'
				Lambda=Set.lambdaS1*Cell.ExternalLambda;
			elseif face.InterfaceType == 'Cell-Cell'
				Lambda=Set.lambdaS2*Cell.InternalLambda;
			elseif face.InterfaceType == 'Bottom'
				Lambda=Set.lambdaS3*Cell.SubstrateLambda;
            else
                disp('ERROR');
			end
			fact0=fact0+Lambda*face.Area;
		end
		fact=fact0/Cell.Area0^2;
        for f=1:length(Cell.Faces)
			face = Cell.Faces(f);
			Tris=Cell.Faces(f).Tris;
			if face.InterfaceType == 'Top'
				Lambda=Set.lambdaS1*Cell.ExternalLambda;
			elseif face.InterfaceType == 'Cell-Cell'
				Lambda=Set.lambdaS2*Cell.InternalLambda;
			elseif face.InterfaceType == 'Bottom'
				Lambda=Set.lambdaS3*Cell.SubstrateLambda;
			end
            for t = 1:length(Tris)
				y1 = Ys(Tris(t).Edge(1),:);
				y2 = Ys(Tris(t).Edge(2),:);
                if length(Tris) == 3
					y3 = Ys(Tris(t+1).Edge(2),:);
					n3 = Cell.globalIds(Tris(t+1).Edge(2));
				else
					y3 = Cell.Faces(f).Centre;
					n3 = Cell.Faces(f).globalIds;
                end
				nY = [Cell.globalIds(Tris(t).Edge)', n3];
				if Geo.Remodelling
					if ~any(ismember(nY,Geo.AssemblegIds))
                        if length(Tris) == 3
                            break
                        else
                		    continue
                        end
					end
				end
				[gs,Ks,Kss]=gKSArea(y1,y2,y3);
				gs=Lambda*gs;
            	ge=Assembleg(ge,gs,nY);
				Ks=fact*Lambda*(Ks+Kss);
				K = AssembleK(K,Ks,nY);
				if length(Tris) == 3
					break
				end
            end
        end
		g=g+ge*fact;
		K=K+(ge)*(ge')/(Cell.Area0^2);
    	EnergyS=EnergyS+ (1/2)*fact0*fact;
	end
end

