function [Geo, Geo_n, Geo_0, Dofs] = ablateCells(Geo, Geo_n, Geo_0, Dofs, Set, t, numStep)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TInitAblation <= t
    if isempty(Geo.cellsToAblate) == 0
        Geo.log = sprintf('%s ---- Performing ablation\n', Geo.log);
        uniqueDebrisCell = Geo.cellsToAblate(1);
        Geo.Cells(uniqueDebrisCell).AliveStatus = 0;
        Geo.Cells(uniqueDebrisCell).ExternalLambda = Set.lambdaSFactor_Debris;
        Geo.Cells(uniqueDebrisCell).InternalLambda = Set.lambdaSFactor_Debris;
        Geo.Cells(uniqueDebrisCell).SubstrateLambda = Set.lambdaSFactor_Debris;
        
        Geo.Cells(uniqueDebrisCell).X = mean(vertcat(Geo.Cells(Geo.cellsToAblate).X));
        
        remainingDebrisCells = setdiff(Geo.cellsToAblate, uniqueDebrisCell);
        for debrisCell = remainingDebrisCells
            [Geo] = CombineTwoNodes(Geo, Set, [uniqueDebrisCell debrisCell]);
            [Geo_n] = CombineTwoNodes(Geo_n, Set, [uniqueDebrisCell debrisCell]);
            [Geo_0] = CombineTwoNodes(Geo_0, Set, [uniqueDebrisCell debrisCell]);
        end
        
        Geo   = Rebuild(Geo, Set);
        Geo   = BuildGlobalIds(Geo);
        Geo   = UpdateMeasures(Geo);
        
        Geo_n = Rebuild(Geo_n, Set);
        Geo_n = BuildGlobalIds(Geo_n);
        Geo_n = UpdateMeasures(Geo_n);
        
        Geo_0 = Rebuild(Geo_0, Set);
        Geo_0 = BuildGlobalIds(Geo_0);
        
        Geo.cellsToAblate = [];
        
        if Set.Substrate == 1
            Dofs = GetDOFsSubstrate(Geo, Set);
        else
            Dofs = GetDOFs(Geo, Set);
        end
        
        PostProcessingVTK(Geo, Geo_0, Set, numStep)
    end
end
end
