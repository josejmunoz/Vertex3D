function [Geo] = ablateCells(Geo, Set, t)
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
        
        remainingDebrisCells = setdiff(Geo.cellsToAblate, uniqueDebrisCell);
        for debrisCell = remainingDebrisCells
            [Geo, Tnew, Ynew] = CombineTwoNodes(Geo, Set, [uniqueDebrisCell debrisCell], [], [])
            Geo.Cells(debrisCell).AliveStatus = [];
            
        end
        
        Geo   = Rebuild(Geo, Set);
        Geo   = BuildGlobalIds(Geo);
        Geo   = UpdateMeasures(Geo);
        
        Geo.cellsToAblate = [];
    end
end
end
