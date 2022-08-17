Set.InputGeo = 'Voronoi';

Set.TotalCells = 20;
Set.CellHeight = 0.2;

Set.Ablation = 1;
Geo.cellsToAblate = [1];
 
Set.lambdaV = 15;

Set.InPlaneElasticity = true;
Set.mu_bulk	= 3000;
Set.lambda_bulk	= 2000; 

Set.tend=42; % 72 = 70 minutes (60 after ablation)
Set.Nincr=Set.tend*80;
Set.TInitAblation = 0.001;
Set.TEndAblation = 2;
 
Set.Substrate  = 2;
Set.kSubstrate = 0.0001;
Set.SubstrateZ = -1;

Set.Contractility = 1;
Set.cLineTension = 0.1;

Set.Remodelling = 1;
Set.RemodelingFrequency = 0.05;
Set.lambdaB = 1e-10;
Set.RemodelTol = Set.lambdaB;

Set.BC = 2;
Set.TStartBC = Set.tend;
Set.TStopBC = Set.tend;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -100;

Set.lambdaS1 = 0.1;
Set.lambdaS2 = 0.000001;
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling_Voronoi';