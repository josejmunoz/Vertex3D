function [Set, CellInput] = updateParametersOnTime(currentT, Set, Cell, CellInput)
%UPDATEPARAMETERSONTIME Summary of this function goes here
%   Detailed explanation goes here

if any(Cell.DebrisCells)
    %% Contractility time dependent
    if currentT > (Set.TAblation + Set.initMidEndContractilityTime_PurseString(2))
        Set.cPurseString = Set.cPurseString_MidEndTimeDependent(currentT);
    else
        Set.cPurseString = Set.cPurseString_InitMidTimeDependent(currentT);
    end
    
    if currentT > (Set.TAblation + Set.initMidEndContractilityTime_LateralCables(2))
        Set.cLateralCables = Set.cLateralCables_MidEndTimeDependent(currentT);
    else
        Set.cLateralCables = Set.cLateralCables_InitMidTimeDependent(currentT);
    end
    
    %% Volume degradation & Surface Cell-Debris/Debris
    if currentT < (Set.TAblation + Set.TToCompleteAblation)
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

