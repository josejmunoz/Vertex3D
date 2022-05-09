function [Geo, Set] = ablateCells(Geo, Set, t)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TInitAblation <= t
    if isempty(Set.cellsToAblate) == 0
        disp('---- Performing ablation');
        debrisCells = ismember(find(cellfun(@isempty, {Geo.Cells.AliveStatus})==0), Set.cellsToAblate); %Improve!!!
        Geo.Cells(debrisCells).debris = true;
        Geo.Cells(debrisCells).ExternalLambda = Set.lambdaSFactor_Debris;
        Geo.Cells(debrisCells).InternalLambda = Set.lambdaSFactor_Debris;
        Geo.Cells(debrisCells).SubstrateLambda = Set.lambdaSFactor_Debris;
        
        Set.cellsToAblate = [];
    end
end
end
