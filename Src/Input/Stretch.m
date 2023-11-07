Geo = struct();
Geo.nx = 5;
Geo.ny = 5;
Geo.nz = 1;
Set.tend=300;
Set.Nincr=300;
Set.BC = 1;
Set.dx = 2; % compression only (2 for stretching)

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.8; % compression only. 0.8 for stretch
Set.VPrescribed = 1.5;
Set.VFixd = -1.5;
Set.ApplyBC=true;

Set.InputGeo = 'Bubbles';
Set.VTK = false;

Set.OutputFolder='Result/Stretch';