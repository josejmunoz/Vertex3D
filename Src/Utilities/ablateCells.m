function [Geo, Set] = ablateCells(Geo, Set, t)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TInitAblation <= t
    if isempty(Set.cellsToAblate) == 0
        disp('---- Performing ablation');
        Cell = Cell.AblateCells(Set.cellsToAblate);
        Set.cellsToAblate = [];
        Geo.Cells(DebrisCells).ExternalLambda = Set.LambdaSFactor_Debris;
        Geo.Cells(DebrisCells).InternalLambda = Set.LambdaSFactor_Debris;
        Geo.Cells(DebrisCells).SubstrateLambda = Set.LambdaSFactor_Debris;
    end
end
end
