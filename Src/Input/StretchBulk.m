Geo.nx = 1;
Geo.ny = 3;
Set.lambdaV = 5;

Set.tend=300;
Set.Nincr=300;
Set.BC = 1;
Set.dx = 2; 

Set.mu_bulk     = 0; % deformation term
Set.lambda_bulk = 0.01; % volume term
Set.InPlaneElasticity = true;

Set.lambdaS1 = 1;
Set.lambdaS2 = 0.8; 
Set.OutputFolder='Result/StretchBulk';