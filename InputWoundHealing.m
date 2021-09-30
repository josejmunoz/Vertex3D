%% geometry
Set.InputSegmentedImage = 'InputImage_dWP3.bmp';

Set.CellHeight = 15; %Microns
Set.zScale = 19.23; %MicronsXY-MicronsZ relation
Set.EllipseFitDiameter = 1; %Microns of a fitted ellipsed in Rob's Wing Discs
Set.AvgCellArea = pi * (Set.EllipseFitDiameter/2)^2; %Microns
Set.CellHeight = (Set.CellHeight * Set.zScale) / Set.AvgCellArea;
Set.TotalCells = 225; %Aim 225
%Set.TotalCells = 40;

%Set.e=4;  % Example Number look in Geo\Example.m 
Set.Method=1;
% Tuning parameters
Set.s=1.5;
Set.f=Set.s/2;

%%  Mechanics
%---------- Volume
Set.lambdaV=20;
Set.lambdaV_Debris=0.01;

%---------- Surface
% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
Set.SurfaceType=4;
% external 
Set.lambdaS1=10;
% Cell-Cell 
Set.lambdaS2=1;
% Cell-substrate
Set.lambdaS3=Set.lambdaS1;

%---------- Line tension
Set.cLineTension = 1;

%---------- In plane elasticity
Set.InPlaneElasticity = true;
Set.mu_bulk = 3000; % Deformation restriction
Set.lambda_bulk = 2000; %Volume restriction

%--------- Bending 
Set.Bending=false;
%------- Viscosity
Set.nu=0.5;   % this is eta 

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
Set.Substrate = true;
Set.kSubstrate = 50;

%% time
Set.tend=1;
Set.Nincr=1000;


%% Remodeling
Set.Remodelling=false;
Set.RemodelTol=.5e-2;
Set.RemodelingFrequency=Set.tend/Set.Nincr;

%---------- EnergyBarrier
Set.EnergyBarrier=true;
Set.lambdaB=1;
Set.Beta=1;
Set.BarrierTri0 = 5e-2; % CARE!! THIS IS OVERRIDE WITHIN THE CODE: INPUTIMAGE.M
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   

%% Ablating cells
Set.Ablation = true;
%Set.cellsToAblate = findCentralCells(Example(Set.e), 1);
% Aim: Set.cellsToAblate = 1:15;
Set.cellsToAblate = 1:15;
Set.TInitAblation = 0.01;
Set.TEndAblation = 0.05;

%% Contractility
% 0: No contractility
% 1: Lateral cables end-to-end
% 2: Lateral surface contractility
Set.Contractility = 1; 

Set.cPurseString = 3;
Set.Contractility_Variability_PurseString = ([1 1 2.5 2] - 1) * Set.cPurseString;
Set.Contractility_TimeVariability_PurseString = [0 7 16 60]/60*(Set.TEndAblation - Set.TInitAblation);

Set.cLateralCables = 0.5;
Set.Contractility_Variability_LateralCables = ([0.5 1.4 1.4] - 0.5) * Set.cLateralCables;
Set.Contractility_TimeVariability_LateralCables = [0 16 60]/60*(Set.TEndAblation - Set.TInitAblation);

%% Execution parameters
Set.OutputFolder = strcat('Result/cLineTension_', num2str(Set.cLineTension),'_typeOfContractility_', num2str(Set.Contractility),'_cPurseString_', num2str(Set.cPurseString), '_cLateralCables_', num2str(Set.cLateralCables), '_lambdaV_', num2str(Set.lambdaV), '_lambdaS1_', num2str(Set.lambdaS1),'_lambda_S2_', num2str(Set.lambdaS2), '_KSubstrate_', num2str(Set.kSubstrate),'_Remodelling_', num2str(Set.Remodelling),'_confinedXYZ_OuterVertices_NCells_', num2str(Set.TotalCells), '_viscosity_', num2str(Set.nu), '_elasticity_mu_', num2str(Set.mu_bulk), '_elasticity_lambda_', num2str(Set.lambda_bulk));
Set.diary = true;
Set.MaxIter = 400;
Set.tol=1e-10;
Set.Parallel = false;
Set.Sparse = 1; %0: No sparse
                    %1: Sparse matlab
                    %2: Sparse manual
