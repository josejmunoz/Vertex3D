Set.OutputFolder = 'Result/Test' ;%'Result/Ablation_Contractility_0.01_NoRemodel_S4_0.5_3x3';
Set.diary = true;
Set.MaxIter = 400;
Set.tol=1e-8;

%% geometry
Set.e=25;  % Example Number look in Geo\Example.m 
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
Set.lambdaS1=1;
% Cell-Cell 
Set.lambdaS2=0.5;
% Cell-substrate
Set.lambdaS3=Set.lambdaS2;
% Cell-GhostCell
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
Set.Remodelling=false;
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=1;

%% time
Set.tend=300;
    
%% Contractility
Set.Contractility = 0;
Set.cContractility = 0.01;
%Set.initEndContractility = [0.001 0.1];
%Set.timeToReachFullContractility = 5;



%% Ablating cells
Set.Ablation = true;
Set.cellsToAblate = findCentralCells(Example(Set.e), 1);
Set.TAblation = 10;
Set.TToCompleteAblation = 100;

%% Contractility
Set.Contractility = 0;

Set.cPurseString = 0.1;
Set.initMidEndContractility_PurseString = ([1 2.5 2] - 1) * Set.cPurseString;
Set.initMidEndContractilityTime_PurseString = [0 18 60]/60*(Set.tend - Set.TAblation);

Set.cLateralCables = 0.05;
Set.initMidEndContractility_LateralCables = ([0.5 1.4 1.4] - 0.5) * Set.cLateralCables;
Set.initMidEndContractilityTime_LateralCables = [0 18 60]/60*(Set.tend - Set.TAblation);
