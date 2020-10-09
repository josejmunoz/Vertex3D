function [Set]=SetDefault(Set)
%% geometry
% Examples of cell centers  
if ~isfield(Set,'e')
    Set.e=1;
end 
% Seeding method==1 % Bounding box   
%         method==2 % DistanceFunction  
if ~isfield(Set,'SeedingMethod')
    Set.SeedingMethod=1;
end 
% Tuning parameters
if ~isfield(Set,'s')
    Set.s=1.5;
end 
if ~isfield(Set,'f')
    Set.f=Set.s/2;
end 
%% Add Substrate
if ~isfield(Set,'Substrate')
    Set.Substrate=false;
    Set.SubstrateZ=0;
end 


%% Time 
Set.tend=180;
Set.Nincr=200;


%%  Mechanics
%---------- Volume
if ~isfield(Set,'lambdaV')
    Set.lambdaV=1;
end 
%---------- Surface
% Set.SurfaceType=1 : Surface-Energy based on the whole Cell-area 
% -         Set.A0eq0=true\false; 
% -         Set.lambdaS=1>0;

% Set.SurfaceType=2 : Surface-Energy based on the Face-area  
% -         Set.lambdaS=1>0;

% Set.SurfaceType=3 : Surface-Energy based on the Triangle-area 
% -         Set.lambdaS=1>0;

% Set.SurfaceType=4 : Surface-Energy based on the whole cell area differential adhsion
%         % - external         Set.lambdaS1=1>0; 
%                                 - Set.LambdaS1CellFactor(:,2)=ones(Celln,1);
%         % - Cell-Cell        Set.lambdaS2=.5>0;
%                                 - Set.LambdaS2CellFactor(:,2)=ones(Celln,1);
%         % - Cell-substrate   Set.lambdaS3=.5>0;
%                                  - Set.LambdaS3CellFactor(:,2)=ones(Celln,1);
if ~isfield(Set,'SurfaceType')
    Set.SurfaceType=1;
    Set.A0eq0=true; 
end 
if ~isfield(Set,'LambdaS1CellFactor')
    Set.LambdaS1CellFactor=[];
end 
if ~isfield(Set,'LambdaS2CellFactor')
    Set.LambdaS2CellFactor=[];
end 
if ~isfield(Set,'LambdaS3CellFactor')
    Set.LambdaS3CellFactor=[];
end 

%---------- EnergyBarrier
% Set.EnergyBarrier=true;
% Set.lambdaB=5;
% Set.Beta=1;  
% WBexp =exp( lambdaB*  ( 1 - Set.Beta*At/At0 )  );   
%  At0 is hard coded as (At0=1e-3) so,  WBexp =exp( lambdaB*  ( 1 - 1000*Set.Beta*At)  );
if ~isfield(Set,'EnergyBarrier')
   Set.EnergyBarrier=true;
end 
if ~isfield(Set,'lambdaB')
     Set.lambdaB=5;
end 
if ~isfield(Set,'Beta')
   Set.Beta=true;
end 

%--------- Bending
if ~isfield(Set,'Bending')    
   Set.Bending=false;
end
if ~isfield(Set,'lambdaBend')
     Set.lambdaBend=0.01;
end
if ~isfield(Set,'BendingAreaDependent')
    Set.BendingAreaDependent=true;
end 

%------- Viscosity
if ~isfield(Set,'nu')
    Set.nu=0.05;   % this is eta 
end 

%% Remodling 
if ~isfield(Set,'RemodelTol')  
    Set.RemodelTol=.5e-6;
end 
if ~isfield(Set,'RemodelingFrequency')  
    Set.RemodelingFrequency=2;
end 

%% Solution 
% ------- Tolerance
Set.tol=1e-10;
% ------- Maximum iteration
Set.MaxIter=30;
Set.Parallel=false;

%% Boundary Condition and loading setting 

% Set.BC=1;  %  Stretch
%     -Set.VFixd=-.5;
%     -Set.VPrescribed=2.5;
%     -Set.dx=1;
%     -Set.TStratBC=20;  %30  
%     -Set.TStopBC=100;
%     
% Set.BC=2;  %  Compression
%        -Set.VFixd=0;
%        -Set.dx=1;
%        -Set.TStratBC=20;  %30  
%        -Set.TStopBC=100;

% Set.BC =~ 1,2;  % Substrate



if ~isfield(Set,'BC') && ~Set.Substrate
    Set.BC=1;
    Set.VFixd=-.5;
    Set.VPrescribed=2.5;
    Set.dx=2;
    Set.TStartBC=20;  %30  
    Set.TStopBC=100;
elseif  Set.Substrate
    Set.BC=nan;
end 

end
       

