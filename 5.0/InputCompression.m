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

%---------- EnergyBarrier
Set.EnergyBarrier=true;
Set.lambdaB=5;
Set.Beta=1;  
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   

%--------- Bending 
Set.Bending=false;
%------- Viscosity
Set.nu=0.05;   % this is eta 

%% Remodeling 
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=2;

%% time
tend=20;
Set.Nincr=10;

%%  Boundary Displacement 
Set.BC=2;  %  Compression
    Set.VFixd=0;
    Set.dx=1;
    Set.TStartBC=20;  %30  
    Set.TStopBC=200;
    
%% Contractility
Set.Contractility = true;
Set.cContractility = 0.3;



%%
    Set.LocalViscosityEdgeBased=true; 
    Set.nu_Local_EdgeBased=0.0001;
    Set.LocalViscosityOption=2;



