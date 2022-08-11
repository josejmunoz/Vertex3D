Set.InputGeo = 'Voronoi';

Set.TotalCells = 20;

Set.Ablation = 1;
Geo.cellsToAblate = [1:3];
 
Set.lambdaV = 10;

% Set.InPlaneElasticity = true;
% Set.mu_bulk	= 300;
% Set.lambda_bulk	= 200; 

Set.tend=42; % 72 = 70 minutes (60 after ablation)
%Set.Nincr=Set.tend*100000;
Set.Nincr=Set.tend*10;
Set.TInitAblation = 0.5;
 
% Set.Substrate  = 2;
% Set.kSubstrate = 0.01;
% Set.SubstrateZ = -0.9;

Set.Contractility = 0;
Set.cLineTension = 0.001;

Set.Remodelling = 1;
Set.RemodelingFrequency = 0.05;
Set.lambdaB = 0.0000001;
Set.RemodelTol = 1e-07;

Set.BC = 2;
Set.TStartBC = Set.tend;
Set.TStopBC = Set.tend;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -1;

Set.lambdaS1 = 0.1;
Set.lambdaS2 = 0.001;
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling';