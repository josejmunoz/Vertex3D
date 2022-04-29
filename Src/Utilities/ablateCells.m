function [Geo, Set] = ablateCells(Geo, Set, t)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TInitAblation <= t
    if isempty(Set.cellsToAblate) == 0
        disp('---- Performing ablation');
        debrisCells = ismember(obj.Int, Set.cellsToAblate);
        Geo.Cells(debrisCells).debris = true;
        Geo.Cells(debrisCells).ExternalLambda = Set.LambdaSFactor_Debris;
        Geo.Cells(debrisCells).InternalLambda = Set.LambdaSFactor_Debris;
        Geo.Cells(debrisCells).SubstrateLambda = Set.LambdaSFactor_Debris;
        
        Set.cellsToAblate = [];
    end
end
end
