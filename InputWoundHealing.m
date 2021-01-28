Set.OutputFolder = 'Result/CompressionAblationNoContractilityRemodel_3x3';
Set.diary = true;

%% geometry
Set.e=4;  % Example Number look in Geo\Example.m 
Set.Method=1;
% Tuning parameters
Set.s=1.5;
Set.f=Set.s/2;

%%  Mechanics
%---------- Volume
Set.lambdaV=5;

%---------- Surface
% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
Set.SurfaceType=4;
% external 
Set.lambdaS1=1;
% Cell-Cell 
Set.lambdaS2=.5;
% Cell-substrate
Set.lambdaS3=.5;
% Cell-GhostCell
Set.lambdaS4=1;

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

%% Remodeling
Set.Remodelling=true;
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=2;

%% time
Set.tend=300;
Set.Nincr=300;
    
%% Contractility
Set.Contractility = 0;
Set.cContractility = 0;

%% Ablating cells
Set.Ablation = true;
Set.cellsToAblate = findCentralCells(Example(Set.e), 1);




