Set.InputGeo = 'Bubbles';

Geo.nx = 1; 
Geo.ny = 3; 
Geo.nz = 1;

%Set.TotalCells = 40;

Set.Ablation = 0;
Geo.cellsToAblate = [5];
 
Set.lambdaV = 20; 
 
Set.tend=200; 
Set.Nincr=Set.tend*40; 
 
Set.Substrate = 0;

Set.Contractility = false;
Set.cLineTension = 0;

Set.Remodelling = 1;
Set.RemodelingFrequency = 0.05;
Set.lambdaB = 0.001;
Set.RemodelTol = Set.lambdaB*0.5;

Set.BC = 2;
Set.TStartBC = 50;
Set.TStopBC = 10;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -1;

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.1; % compression only. 0.8 for stretch
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling_Bubbles';