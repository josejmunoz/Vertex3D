Geo = struct();
Geo.nx = 10;
Geo.ny = 10;
Set.tend=300;
Set.Nincr=300;
Set.BC = 1;
Set.dx = 2; % compression only (2 for stretching)

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.8; % compression only. 0.8 for stretch
Set.VPrescribed = 1.5;
Set.VFixd = -1.5;
Set.ApplyBC=true;

Set.OutputFolder='Result/Stretch';