function [CellInput, Set]=InitializeInput(Cell,Set)
% Initialize Lmabda-reduction-factor for each cell-inerface (check KgSurfaceCellBasedAdhesionParallel.m)
% Initialize propulsion force for each cell (check gPropulsion.m)
% Initialize contractility in time

if Set.SurfaceType==4
    aux1=Set.LambdaS1CellFactor;
    aux2=Set.LambdaS2CellFactor;
    aux3=Set.LambdaS3CellFactor;
    aux4=Set.LambdaS4CellFactor; 
    CellInput.LambdaS1Factor=ones(Cell.n,1);
    CellInput.LambdaS2Factor=ones(Cell.n,1);
    CellInput.LambdaS3Factor=ones(Cell.n,1);
    CellInput.LambdaS4Factor=ones(Cell.n,1);
    if ~isempty(aux1)
        for i=1:size(aux1,1)
            CellInput.LambdaS1Factor(aux1(i,1))=aux1(i,2);
        end
    end
    
    if ~isempty(aux2)
        for i=1:size(aux2,1)
            CellInput.LambdaS2Factor(aux2(i,1))=aux2(i,2);
        end
    end
    
    if ~isempty(aux3)
        for i=1:size(aux3,1)
            CellInput.LambdaS3Factor(aux3(i,1))=aux3(i,2);
        end
    end
    
    if ~isempty(aux4)
        for i=1:size(aux4,1)
            CellInput.LambdaS4Factor(aux4(i,1))=aux4(i,2);
        end
    end   
else 
    CellInput.LambdaS1Factor=[];
    CellInput.LambdaS2Factor=[];
    CellInput.LambdaS3Factor=[];
    CellInput.LambdaS4Factor=[];
end 

if Set.Propulsion
    CellInput.Propulsion=zeros(Cell.n,3);
    CellInput.Propulsion(:,1:2)=(rand(Cell.n,2)-0.5)./50;
else 
    CellInput.Propulsion=[];
end 

if isempty(Set.initMidEndContractility_PurseString) == 0
    initT_Contractility = Set.initMidEndContractilityTime_PurseString(1);
    midT_Contractility = Set.initMidEndContractilityTime_PurseString(2);
    endT_Contractility = Set.initMidEndContractilityTime_PurseString(3);
    
    initC_Contractility = Set.initMidEndContractility_PurseString(1);
    midC_Contractility = Set.initMidEndContractility_PurseString(2);
    endC_Contractility = Set.initMidEndContractility_PurseString(3);
    
    [Set.cPurseString_InitMidTimeDependent] = functionVariableOnTime(initC_Contractility, midC_Contractility, initT_Contractility, midT_Contractility, Set);
    
    [Set.cPurseString_MidEndTimeDependent] = functionVariableOnTime(midC_Contractility, endC_Contractility, midT_Contractility, endT_Contractility, Set);
end

if isempty(Set.initMidEndContractility_LateralCables) == 0
    initT_Contractility = Set.initMidEndContractilityTime_LateralCables(1);
    midT_Contractility = Set.initMidEndContractilityTime_LateralCables(2);
    endT_Contractility = Set.initMidEndContractilityTime_LateralCables(3);
    
    initC_Contractility = Set.initMidEndContractility_LateralCables(1);
    midC_Contractility = Set.initMidEndContractility_LateralCables(2);
    endC_Contractility = Set.initMidEndContractility_LateralCables(3);
    
    [Set.cLateralCables_InitMidTimeDependent] = functionVariableOnTime(initC_Contractility, midC_Contractility, initT_Contractility, midT_Contractility, Set);
    
    [Set.cLateralCables_MidEndTimeDependent] = functionVariableOnTime(midC_Contractility, endC_Contractility, midT_Contractility, endT_Contractility, Set);
end

if Set.TToCompleteAblation > 0
    Set.lambdaV_DebrisTime = functionVariableOnTime(Set.lambdaV, Set.lambdaV_Debris, 0, Set.TToCompleteAblation, Set);
    Set.lambdaS4_Time = functionVariableOnTime(Set.lambdaS2, Set.lambdaS4, 0, Set.TToCompleteAblation, Set);
end


end 
