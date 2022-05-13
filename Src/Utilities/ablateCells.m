function [Geo] = ablateCells(Geo, Set, t)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TInitAblation <= t
    if isempty(Geo.cellsToAblate) == 0
        disp('---- Performing ablation');
        for debrisCell = Geo.cellsToAblate
            Geo.Cells(debrisCell).AliveStatus = 0;
            Geo.Cells(debrisCell).ExternalLambda = Set.lambdaSFactor_Debris;
            Geo.Cells(debrisCell).InternalLambda = Set.lambdaSFactor_Debris;
            Geo.Cells(debrisCell).SubstrateLambda = Set.lambdaSFactor_Debris;
        end
        Geo.cellsToAblate = [];
    end
end
end
