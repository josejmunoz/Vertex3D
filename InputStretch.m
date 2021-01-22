% Meaning of parameters is in SetDefault.m
%
%% geometry
Set.e=1;  % Example Number look in Geo\Example.m 
Set.Method=1; % Method for seeding the boundary nodes 
% Tuning parameters
Set.s=1.5; % Average cell size 
Set.f=Set.s/2; % Average distance from vertex to cell centre

%%  Mechanics
%---------- Volume
Set.lambdaV=5; % Penalisation value for volume

%---------- Surface
% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhesion
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
%------- Global Viscosity coefficient 
Set.nu=0.05;   % this is eta 

%% Remodeling 
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=2;

%% time
tend=300;
Set.Nincr=300;

%%  Boundary Displacement 
Set.BC=1;  %  Stretching
    Set.VFixd=-1.5;      % Vertices with y coordinate < VFixs will not move
    Set.VPrescribed=1.5; % Vertices with y coordinate > VPrescribed will have displacement imposed
    Set.dx=2;         % Incremetnal stretching displcacement on side
    Set.TStartBC=20;  % Time where stretch is initiated
    Set.TStopBC=200;  % Time where stretch is stopped



% %% time Example 1 Stratch 
% tend=200;
% Set.Tstretch=20;  %30  
% Set.Trelax=200;
% Set.Nincr=300;
% 
% 
% %%  Boundary Displacement 
% Set.dx=2;

