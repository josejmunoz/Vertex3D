%% geometry
Set.e=7;  % Example Number look in Geo\Example.m 
Set.Method=1;
% Tuning parameters
Set.s=1.5;
Set.f=Set.s/2;

%%  Mechanics
%% Volume
Set.lambdaV=20;

%% Surface acinar cells
% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
Set.SurfaceType=4;
% external 
Set.lambdaS1=2;
% Cell-Cell 
Set.lambdaS2=1;
% Cell-substrate
Set.lambdaS3=2;

%Surface acinar cell to create inequality
Set.LambdaS1CellFactor=[8 2];
Set.LambdaS2CellFactor=[8 2];
Set.LambdaS3CellFactor=[8 2];

%% Surface ductal cells
Set.LambdaS1CellFactor=[10 2];
Set.LambdaS2CellFactor=[10 2];
Set.LambdaS3CellFactor=[10 2];

%% Substrate
Set.Substrate = false;
% Set.kSubstrate = 0.01;
% Set.z0Substrate = 3;

%% Cell movement (ductal cell)
Set.CellMovement = true;
Set.MovementStrength = 0.01;
Set.DestinationPoint = [0, 0, -10];

%% EnergyBarrier
Set.EnergyBarrier=true;
Set.lambdaB=5;
Set.Beta=1;  
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   

%% Bending 
Set.Bending=false;

%% Viscosity
Set.nu=0.1;   % this is eta

%% Contractility
Set.Contractility = 0; 
Set.cLineTensionApical = 0.0001;
Set.cLineTensionBasal = 0.0001;
Set.cLineTensionLateral = 0.0001;

%% time
Set.tend=10;
Set.Nincr=Set.tend*20;

%% Remodeling
Set.Remodelling=true;
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=Set.tend / Set.Nincr;

%%  Boundary Displacement 
Set.BC=2;  %  Compression
    Set.VFixd=0;
    Set.dx=1;
    Set.TStartBC=20;  %30  
    Set.TStopBC=200;

Set.cellsToAblate = 0;
    
%%
Set.batchProcessing = 0;




