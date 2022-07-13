Geo.nx = 1; 
Geo.ny = 3; 
Geo.nz = 1;

Set.Ablation = 0;
Geo.cellsToAblate = [5];
 
Set.lambdaV = 10; 
 
Set.tend=200; 
Set.Nincr=400; 
 
Set.Substrate  = 2;
Set.kSubstrate = 0.01;
Set.SubstrateZ = -0.9;

Set.Contractility = true;
Set.cLineTension = 0.001;

%Set.Remodelling = 0;
Set.RemodelingFrequency = 0.5;

Set.BC = 1;
Set.dx = 2; % compression only (2 for stretching)

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.8; % compression only. 0.8 for stretch
Set.VPrescribed = 1.5;
Set.VFixd = -1.5;
Set.ApplyBC=true;
 
Set.OutputFolder='Result/Remodelling';