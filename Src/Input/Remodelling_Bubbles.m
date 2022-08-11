Set.InputGeo = 'Bubbles';

Geo.nx = 3; 
Geo.ny = 1; 
Geo.nz = 1;

%Set.TotalCells = 40;

Set.Ablation = 0;
Geo.cellsToAblate = [5];
 
Set.lambdaV = 10; 
 
Set.tend=200; 
Set.Nincr=800; 
 
Set.Substrate = 0;

Set.Contractility = true;
Set.cLineTension = 0.001;

Set.Remodelling = 1;
Set.RemodelingFrequency = 0.05;
Set.lambdaB = 0.0000001;
Set.RemodelTol = 6e-08;

Set.BC = 2;
Set.TStartBC = Set.tend;
Set.TStopBC = Set.tend;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -1;

Set.lambdaS1 = 0.1;
Set.lambdaS2 = 0.01; % compression only. 0.8 for stretch
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling';