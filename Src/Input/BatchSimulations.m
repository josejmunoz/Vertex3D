Set.InputGeo = 'Voronoi';

Set.TotalCells = 100;
Set.CellHeight = 15;

Geo.cellsToAblate = [1:5];
 
Set.lambdaV = 20;

Set.InPlaneElasticity = true;
Set.mu_bulk	= 1000;
Set.lambda_bulk	= 500;

Set.tend=61; % 72 = 70 minutes (60 after ablation)
Set.Nincr=Set.tend*2;
Set.TInitAblation = 1;
Set.TEndAblation = Set.tend;
 
Set.Substrate  = 2;
Set.kSubstrate = 500;
Set.SubstrateZ = -Set.CellHeight/2;

Set.Contractility = 1;
Set.cLineTension = 0.03;

Set.Remodelling = 1;
Set.RemodelingFrequency = (Set.tend/Set.Nincr)*2;
Set.lambdaB = 1e-4;
Set.RemodelTol = 0;
Set.RemodelStiffness = 0.9;

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

