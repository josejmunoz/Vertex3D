function [Set, CellInput] = updateParametersOnTime(currentT, Set, Cell, CellInput)
%UPDATEPARAMETERSONTIME Summary of this function goes here
%   Detailed explanation goes here

if any(Cell.DebrisCells)
    %% Contractility time dependent
    % Purse string
    if isempty(Set.Contractility_TimeVariability_PurseString) == 0
        achievedTimes = find(currentT >= (Set.TInitAblation + Set.Contractility_TimeVariability_PurseString));
        Set.cPurseString = Set.cPurseString_TimeDependent{achievedTimes(end)}(currentT);
    end
    
    % Lateral cables
    if isempty(Set.Contractility_TimeVariability_LateralCables) == 0
        achievedTimes = find(currentT >= (Set.TInitAblation + Set.Contractility_TimeVariability_LateralCables));
        Set.cLateralCables = Set.cLateralCables_TimeDependent{achievedTimes(end)}(currentT);
    end

    %% Volume degradation & Surface Cell-Debris/Debris
    if currentT < (Set.TInitAblation + Set.TEndAblation)
        Set.lambdaV_Debris = Set.lambdaV_DebrisTime(currentT);
        Set.lambdaS4 = Set.lambdaS4_Time(currentT);
        CellInput.LambdaS1Factor(Cell.DebrisCells) = Set.LambdaS1FactorDebris_Time(currentT);
    end
else
    Set.cPurseString = 0;
    Set.cLateralCables = 0;
    Set.lambdaV_Debris = Set.lambdaV;
    Set.lambdaS4 = Set.lambdaS2;
end
end

