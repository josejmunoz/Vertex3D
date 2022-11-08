Set.InputGeo = 'Voronoi';

Set.TotalCells = 40;
Set.CellHeight = 0.2;

Set.Ablation = 1;
Geo.cellsToAblate = [1:3];
 
Set.lambdaV = 20;

Set.InPlaneElasticity = true;
Set.mu_bulk	= 1000;
Set.lambda_bulk	= 500;

Set.tend=42; % 72 = 70 minutes (60 after ablation)
Set.Nincr=Set.tend*80;
Set.TInitAblation = 0.001;
Set.TEndAblation = 2;
 
Set.Substrate  = 2;
Set.kSubstrate = 500;
Set.SubstrateZ = -Set.CellHeight/2;

Set.Contractility = 1;
Set.cLineTension = 0.01;

Set.Remodelling = 1;
Set.RemodelingFrequency = 0.02;
Set.lambdaB = 1e-4;
Set.RemodelTol = 0.02;
Set.RemodelStiffness = -0.1;

Set.BC = 2;
Set.TStartBC = Set.tend;
Set.TStopBC = Set.tend;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -100;
Set.Reset_PercentageGeo0 = 0.3; 

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.9;
Set.lambdaS3 = 1;
Set.ApplyBC=true;
 
Set.OutputFolder=strcat('Result/Remodelling_Voronoi_', num2str(Set.TotalCells));