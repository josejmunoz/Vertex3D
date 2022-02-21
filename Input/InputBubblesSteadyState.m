%% geometry
Set.e=7;  % Example Number look in Geo\Example.m 
Set.Method=1;
% Tuning parameters
Set.s=1.5;
Set.f=Set.s/2;

%%  Mechanics
%% Volume
Set.lambdaV=5;

%% Surface acinar cells
% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
Set.SurfaceType=4;
% external 
Set.lambdaS1=1;
% Cell-Cell 
Set.lambdaS2=0.5;
% Cell-substrate
Set.lambdaS3=1;

%% Surface ductal cells
Set.LambdaS1CellFactor=[10 5];
Set.LambdaS2CellFactor=[10 5];
Set.LambdaS3CellFactor=[10 5];

%% Substrate pulling
Set.Substrate = true;
Set.kSubstrate = 1;
Set.z0Substrate = 3;

%% EnergyBarrier
Set.EnergyBarrier=true;
Set.lambdaB=5;
Set.Beta=1;  
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   

%% Bending 
Set.Bending=false;

%% Viscosity
Set.nu=0.05;   % this is eta

%% Contractility
Set.Contractility = 1; Set.cLineTension = 0;

%% Remodeling
Set.Remodelling=false;
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=2;

%% time
Set.tend=10;
Set.Nincr=Set.tend*2;

%%  Boundary Displacement 
Set.BC=2;  %  Compression
    Set.VFixd=0;
    Set.dx=1;
    Set.TStartBC=20;  %30  
    Set.TStopBC=200;

Set.cellsToAblate = 0;
    
%%
Set.batchProcessing = 0;




