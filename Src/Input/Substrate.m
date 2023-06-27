Set.InputGeo = 'Bubbles';

Geo.nx = 3; 
Geo.ny = 1; 
Geo.nz = 1;

Set.TotalCells = 40;

Set.Ablation = 0;
Geo.cellsToAblate = [5];
 
Set.lambdaV = 20; 
 
Set.tend=200; 
Set.Nincr=400; 
 
Set.Substrate  = 2;
Set.kSubstrate = 0.01;
Set.SubstrateZ = -0.9;

Set.Contractility = true;
Set.cLineTension = 0.001;

Set.Remodelling = 0;
Set.RemodelingFrequency = 0.5;
Set.lambdaB = 0.0001; % THE BIGGER THE MORE ENERGY

Set.BC = 2;
Set.dx = 1; % compression only (2 for stretching)
Set.VPrescribed = realmax;
Set.VFixd = -1;

Set.lambdaS1 = 1;
Set.lambdaS2 = 1; % compression only. 0.8 for stretch
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling';