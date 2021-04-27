%% geometry
Set.InputSegmentedImage = 'InputImage_dWP3.bmp';

Set.CellHeight = 35; %Microns
Set.zScale = 19.23; %MicronsXY-MicronsZ relation
Set.AvgCellArea = 5; %Microns
Set.CellHeight = (Set.CellHeight * Set.zScale) / Set.AvgCellArea;
Set.TotalCells = 40;

%Set.e=4;  % Example Number look in Geo\Example.m 
Set.Method=1;
% Tuning parameters
Set.s=1.5;
Set.f=Set.s/2;

%%  Mechanics
%---------- Volume
Set.lambdaV=5;
Set.lambdaV_Debris=0.001;

%---------- Surface
% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
Set.SurfaceType=4;
% external 
Set.lambdaS1=0.5;
% Cell-Cell 
Set.lambdaS2=0.1;
% Cell-substrate
Set.lambdaS3=Set.lambdaS2;
% Cell-DebrisCell
Set.lambdaS4=Set.lambdaS2;

%---------- EnergyBarrier
Set.EnergyBarrier=true;
Set.lambdaB=5;
Set.Beta=1;  
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   

%--------- Bending 
Set.Bending=false;
%------- Viscosity
Set.nu=0.05;   % this is eta 

% Set.Confinement=true;
% Set.ConfinementX1=0.5;    
% Set.ConfinementY1=0.5;
% Set.ConfinementX2=-0.5;
% Set.ConfinementY2=-0.5;

%% Compression or stretching
Set.BC=2; % BC=1: Stretching, BC=2: Compression, BC=nan, substrate extrussion
    Set.VFixd=-1.5;
    Set.VPrescribed=1.5;
    Set.dx=0;
    Set.TStartBC=301;  %30  
    Set.TStopBC=302;
    
%% Substrate
Set.Substrate = false;
Set.kSubstrate = 0;

%% Remodeling
Set.Remodelling=true;
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=1;

%% time
Set.tend=20;
Set.Nincr=400;

%% Ablating cells
Set.Ablation = true;
%Set.cellsToAblate = findCentralCells(Example(Set.e), 1);
Set.cellsToAblate = 1:3;
Set.TInitAblation = 1;
Set.TEndAblation = 5;

%% Contractility
Set.Contractility = true;

Set.cPurseString = 0.1;
Set.Contractility_Variability_PurseString = ([1 1 2.5 2] - 1) * Set.cPurseString;
Set.Contractility_TimeVariability_PurseString = [0 7 16 60]/60*(Set.TEndAblation - Set.TInitAblation);

Set.cLateralCables = 0.1;
Set.Contractility_Variability_LateralCables = ([0.5 1.4 1.4] - 0.5) * Set.cLateralCables;
Set.Contractility_TimeVariability_LateralCables = [0 16 60]/60*(Set.TEndAblation - Set.TInitAblation);

%% Execution parameters
Set.OutputFolder = strcat('Result/cellHeight_', num2str(Set.CellHeight),'_cPurseString_', num2str(Set.cPurseString), '_cLateralCables_', num2str(Set.cLateralCables), '_lambdaV_', num2str(Set.lambdaV), '_lambdaS1_', num2str(Set.lambdaS1),'_lambda_S2_', num2str(Set.lambdaS2),'_KSubstrate_', num2str(Set.kSubstrate),'_Remodelling_', num2str(Set.Remodelling),'_confinedXYZ_OuterVertices_NCells_', num2str(Set.TotalCells), '_viscosity_', num2str(Set.nu));
Set.diary = true;
Set.MaxIter = 400;
Set.tol=1e-10;
Set.Parallel = false;
Set.Sparse = false; %0: No sparse
                    %1: Sparse matlab
                    %2: Sparse manual
