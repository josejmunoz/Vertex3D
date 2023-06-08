Set.InputGeo = 'Voronoi';
Geo.cellsToAblate = [1:3];

Set.InPlaneElasticity = false;
Set.mu_bulk	= 0; %% 1000
Set.lambda_bulk	= 0; %% 500

Set.nu	= 5000;
Set.Nincr = 61*10;

Set.lambdaB	= 5;

Set.lambdaV = 10;

Set.kSubstrate = 100;

Set.cLineTension = 0.5;
Set.cLineTensionMembrane = 0;
Set.purseStringStrength = 2.5;

Set.lambdaS1 = 3;
Set.VTK = false;