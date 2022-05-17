function [g,K,EnergyB]=KgTriEnergyBarrier(Geo,Set)
	% The residual g and Jacobian K of  Energy Barrier
	% Energy  WBexp = exp( Set.lambdaB*  ( 1 - Set.Beta*At/Set.BarrierTri0 )  );

	[g, K] = initializeKg(Geo, Set);
	EnergyB = 0;
	for c=1:Geo.nCells
		if Geo.Remodelling
			if ~ismember(c,Geo.AssembleNodes)
        		continue
			end
		end
		Cell = Geo.Cells(c);
		Ys = Cell.Y;
		lambdaB=Set.lambdaB;
		for f = 1:length(Cell.Faces)
            Face = Cell.Faces(f);
			Tris = Cell.Faces(f).Tris;
			for t = 1:length(Tris)
				fact=-((lambdaB*Set.Beta)/Set.BarrierTri0)* ... % * ...
				                    exp(lambdaB*(1-Set.Beta*Face.Tris(t).Area/Set.BarrierTri0));
				fact2=fact*-((lambdaB*Set.Beta)/Set.BarrierTri0);
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
	        	g=Assembleg(g,gs*fact,nY);	
				Ks=(gs)*(gs')*fact2+Ks*fact+Kss*fact;
				K= AssembleK(K,Ks,nY);
				EnergyB=EnergyB+ exp(lambdaB*(1-Set.Beta*Face.Tris(t).Area/Set.BarrierTri0));
%                 fprintf("%.12f %.12f %.12f %.3f %.3f %3f %d %d %d\n", norm(g), norm(K), EnergyB, y3, c, f, t);
				if length(Tris) == 3
% 					sort(gs), sort(Ks)
					break
				end
			end
		end
	end
end
