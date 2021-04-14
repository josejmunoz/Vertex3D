function [Cell, Set, CellInput] = performAblation(Cell, Set, CellInput, t)
%PERFORMABLATION Summary of this function goes here
%   Detailed explanation goes here
if Set.Ablation == true && Set.TInitAblation <= t
    if isempty(Set.cellsToAblate)==0
        Cell = Cell.AblateCells(Set.cellsToAblate);
        Set.cellsToAblate = [];
        CellInput.LambdaS1Factor(Cell.DebrisCells) = 0.00001;
        CellInput.LambdaS2Factor(Cell.DebrisCells) = 0.00001;
        CellInput.LambdaS3Factor(Cell.DebrisCells) = 0.00001;
        
        %% Smaller time-steps
        disp('Updating time-step after ablation');
        Set.Nincr = Set.Nincr*20;
        Set.dt0=Set.tend/Set.Nincr;
        Set.dt=Set.dt0;
    end
end
end

