Set.InputGeo = 'Voronoi';
Geo.cellsToAblate = [1:3];

Set.InPlaneElasticity = false;
Set.mu_bulk	= 0; %% 1000
Set.lambda_bulk	= 0; %% 500

Set.nu	= 5000;
Set.lambdaV = 10;

Set.kSubstrate = 100;

Set.cLineTension = 0.05;
Set.purseStringStrength = 2.5;
Set.RemodelStiffness = 0.65;

Set.lambdaS1 = 0.5;
Set.lambdaS2 = Set.lambdaS1 * 0.9;
Set.lambdaS3 = Set.lambdaS1;