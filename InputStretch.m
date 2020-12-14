%% geometry
Set.e=1;  % Example Number look in Geo\Example.m 
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
Set.lambdaS2=.8;
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
tend=300;
Set.Nincr=300;

%%  Boundary Displacement 
Set.BC=1;  %  Compression
    Set.VFixd=-.5;
    Set.VPrescribed=2.5;
    Set.dx=2;
    Set.TStartBC=20;  %30  
    Set.TStopBC=200;



% %% time Example 1 Stratch 
% tend=200;
% Set.Tstretch=20;  %30  
% Set.Trelax=200;
% Set.Nincr=300;
% 
% 
% %%  Boundary Displacement 
% Set.dx=2;

