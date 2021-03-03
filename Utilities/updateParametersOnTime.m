function [Set] = updateParametersOnTime(currentT, Set, Cell)
%UPDATEPARAMETERSONTIME Summary of this function goes here
%   Detailed explanation goes here

if any(Cell.GhostCells)
    %% Contractility time dependent
    if currentT > (Set.TAblation + Set.initMidEndContractilityTime_PurseString(2))
        Set.cPurseString = Set.cPurseString_MidEndTimeDependent(currentT) + Set.initMidEndContractility_PurseString(2);
    else
        Set.cPurseString = Set.cPurseString_InitMidTimeDependent(currentT) + Set.initMidEndContractility_PurseString(1);
    end
    
    if currentT > (Set.TAblation + Set.initMidEndContractilityTime_LateralCables(2))
        Set.cLateralCables = Set.cLateralCables_MidEndTimeDependent(currentT) + Set.initMidEndContractility_LateralCables(2);
    else
        Set.cLateralCables = Set.cLateralCables_InitMidTimeDependent(currentT) + Set.initMidEndContractility_LateralCables(1);
    end
    
    %% Volume degradation
    if currentT < (Set.TAblation + Set.TToCompleteAblation)
        Set.lambdaV_Debris = Set.lambdaV_DebrisTime(currentT);
    end
else
    Set.cPurseString = 0;
    Set.cLateralCables = 0;
    Set.lambdaV_Debris = Set.lambdaV;
end
end

