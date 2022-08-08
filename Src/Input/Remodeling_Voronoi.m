Set.InputGeo = 'Voronoi';

Set.TotalCells = 10;

Set.Ablation = 0;
Geo.cellsToAblate = [5];
 
Set.lambdaV = 5;

Set.InPlaneElasticity = true;
Set.mu_bulk	= 300;
Set.lambda_bulk	= 200;
 
Set.tend=200; 
Set.Nincr=400; 
 
Set.Substrate  = 2;
Set.kSubstrate = 0.01;
Set.SubstrateZ = -0.9;

Set.Contractility = true;
Set.cLineTension = 0.001;

Set.Remodelling = 0;
Set.RemodelingFrequency = 0.5;
Set.lambdaB = 0.0001;

Set.BC = 2;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -1;

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.001;
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling';