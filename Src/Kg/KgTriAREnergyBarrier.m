function [g,K,EnergyB]=KgTriAREnergyBarrier(Geo,Set)
%%KGTRIARENERGYBARRIER Penalise bad aspect ratios
% The residual g and Jacobian K of  Energy Barrier
% Energy  WB =
    [g, K] = initializeKg(Geo, Set);
    EnergyB = 0;
    
    for c=1:Geo.nCells
        if Geo.Remodelling
            if ~ismember(c,Geo.AssembleNodes)
                continue
            end
        end

        if isempty(Geo.Cells(c).AliveStatus) || Geo.Cells(c).AliveStatus ~= 1
            continue
        end

        Cell = Geo.Cells(c);
        Ys = Cell.Y;
        for f = 1:length(Cell.Faces)
            Face = Cell.Faces(f);
            Tris = Cell.Faces(f).Tris;
            for t = 1:length(Tris)
                n3 = Cell.Faces(f).globalIds;
				nY = [Cell.globalIds(Tris(t).Edge)', n3];
                if Geo.Remodelling
                    if ~any(ismember(currentTet_ids,Geo.AssemblegIds))
                        continue
                    end
                end
                
                y1 = Ys(Tris(t).Edge(1),:);
                y2 = Ys(Tris(t).Edge(2),:);
                y3 = Cell.Faces(f).Centre;
                
                g=Assembleg(g,gs*fact,nY);
                K= AssembleK(K,Ks,nY);
                EnergyB=EnergyB + exp(lambdaB*(1-Set.Beta*Face.Tris(t).Area/Set.BarrierTri0));
            end
        end
    end
end
