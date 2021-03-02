function [Set] = updateContractilityOnTime(currentT, Set, Cell)
%UPDATECONTRACTILITYONTIME Summary of this function goes here
%   Detailed explanation goes here
%% Contractility time dependent
if any(Cell.GhostCells)
    if currentT < (Set.TAblation + Set.initMidEndContractilityTime_PurseString(2))
        Set.cPurseString = Set.cPurseString_MidEndTimeDependent(currentT);
    else
        Set.cPurseString = Set.cPurseString_InitMidTimeDependent(currentT);
    end
else
    Set.cPurseString = 0;
end

if any(Cell.GhostCells)
    if currentT < (Set.TAblation + Set.initMidEndContractilityTime_LateralCables(2))
        Set.cLateralCables = Set.cLateralCables_MidEndTimeDependent(currentT);
    else
        Set.cLateralCables = Set.cLateralCables_InitMidTimeDependent(currentT);
    end
else
    Set.cLateralCables = 0;
end
end

