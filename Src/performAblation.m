function [Cell, Set, CellInput] = performAblation(Cell, Set, CellInput, t)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TAblation <= t
    if isempty(Set.cellsToAblate)==0
        Cell = Cell.AblateCells(Set.cellsToAblate);
        Set.cellsToAblate = [];
        CellInput.LambdaS1Factor(Cell.GhostCells) = 0;
        CellInput.LambdaS2Factor(Cell.GhostCells) = 0;
        CellInput.LambdaS3Factor(Cell.GhostCells) = 0;
    end
    
    if isempty(Set.initEndContractility) == 0
        Set.cContractility = Set.initEndContractility(Set.timeToReachFullContractility);
        Set.timeToReachFullContractility = Set.timeToReachFullContractility - 1;
    end
end
end

