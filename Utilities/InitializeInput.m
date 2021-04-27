function [CellInput, Set]=InitializeInput(Cell, Set, Y)
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

%% Set z0Substrate
Set.z0Substrate = min(Y.DataRow(:,3))*1.01;


%% Contractility
if isempty(Set.Contractility_Variability_PurseString) == 0
    for numTimePoint = 2:length(Set.Contractility_Variability_PurseString)
        [Set.cPurseString_TimeDependent{numTimePoint-1}] = functionVariableOnTime(Set.Contractility_Variability_PurseString(numTimePoint-1), Set.Contractility_Variability_PurseString(numTimePoint), Set.Contractility_TimeVariability_PurseString(numTimePoint-1), Set.Contractility_TimeVariability_PurseString(numTimePoint), Set);
    end
end

if isempty(Set.Contractility_Variability_LateralCables) == 0
    for numTimePoint = 2:length(Set.Contractility_Variability_LateralCables)
        [Set.cLateralCables_TimeDependent{numTimePoint-1}] = functionVariableOnTime(Set.Contractility_Variability_LateralCables(numTimePoint-1), Set.Contractility_Variability_LateralCables(numTimePoint), Set.Contractility_TimeVariability_LateralCables(numTimePoint-1), Set.Contractility_TimeVariability_LateralCables(numTimePoint), Set);
    end
end

if Set.TEndAblation > 0
    Set.lambdaV_DebrisTime = functionVariableOnTime(Set.lambdaV, Set.lambdaV_Debris, Set.TInitAblation, Set.TEndAblation, Set);
    Set.lambdaS4_Time = functionVariableOnTime(Set.lambdaS2, Set.lambdaS4, Set.TInitAblation, Set.TEndAblation, Set);
    %Set.LambdaS1FactorDebris_Time = functionVariableOnTime(0.001, 0.001, 0.001, Set.TEndAblation/10, Set);
end
end 
