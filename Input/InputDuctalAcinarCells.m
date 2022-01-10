%% geometry
Set.e=10;  % Example Number look in Geo\Example.m 
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
% Set.LambdaS1CellFactor=[CellNumber factor];
% Cell-Cell 
Set.lambdaS2=0.2;

%---------- EnergyBarrier
Set.EnergyBarrier=true;
Set.lambdaB=5;
Set.Beta=1;  
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   

%--------- Bending 
Set.Bending=false;
%------- Viscosity
Set.nu=0.2;   % this is eta 

%% Remodeling 
Set.RemodelTol=.5e-6;
Set.RemodelingFrequency=2;

%% time
tend=300;
Set.Nincr=300;
%%  Boundary Displacement 
Set.BC=2;  %  Compression
    Set.VFixd=-4;
    Set.dx=1;
    Set.TStartBC=3;  %30  
    Set.TStopBC=200;


% 
% 
% 
% %% geometry
% Set.e=4;  % Example Number look in Geo\Example.m 
% 
% % Seeding method==1 % Bounding box   
% %         method==2 % DistanceFunction  
% 
% Set.Method=1;
% % Tuning parameters
% Set.s=1.5;
% Set.f=Set.s/2;
% 
% 
% % Set.Method=2;
% % Set.s=2;
% % Set.f=Set.s/1.5;
% 
%     
% 
% 
% %%  Mechanics
% 
% %---------- Volume
% Set.lambdaV=5;
% 
% %---------- Surface
% % Set.SurfaceType=1 : Surface-Energy based on the whole Cell-ar ea  
% % Set.SurfaceType=2 : Surface-Energy based on the Face-area  
% % Set.SurfaceType=3 : Surface-Energy based on the Triangle-area 
% % Set.A0eq0=false;
% Set.A0eq0=true;
% Set.SurfaceType=1;
% Set.lambdaS=.2;
% % Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
% Set.SurfaceType=4;
% % external 
% Set.lambdaS1=1;
% % Cell-Cell 
% Set.lambdaS2=.5;
% % Cell-substrate
% Set.lambdaS3=.5;
% 
% %---------- EnergyBarrier
% Set.EnergyBarrier=true;
% Set.lambdaB=5;
% Set.Beta=1;  
% % WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   
% 
% 
% %--------- Bending 
% Set.lambdaBend=.00;
% Set.BendingAreaDependent=true;
% %------- Viscosity
% Set.nu=0.05;   % this is eta 
% Set.nu0=Set.nu;
% 
% %% Remodeling 
% Set.RemodelTol=.5e-6;
% 
% 
% %% Solution 
% % ------- Tolerance
% Set.tol=1e-10;
% % ------- Maximum iteration
% Set.MaxIter0=30;
% Set.MaxIter=Set.MaxIter0;
% 
% 
% %% time Example 1 Stratch 
% tend=180;
% Set.Nincr=200;
% 
% Set.Substrate=.5;
% 

