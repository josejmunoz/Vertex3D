Set.InputGeo = 'Voronoi';
% 40 cells; 3 cells to ablate
% 110 cells; 5 cells to ablate
Set.TotalCells = 110;
Geo.cellsToAblate = [1:5];

Set.InPlaneElasticity = false;
Set.mu_bulk	= 0; %% 1000
Set.lambda_bulk	= 0; %% 500

Set.nu	= 5000;
Set.Nincr = 61*2;

Set.lambdaB	= 1; % Smaller number, more energy
Set.lambdaR	= 0.1; 

Set.lambdaV = 1000;

Set.kSubstrate = 100;

Set.cLineTension = 0.005;
Set.purseStringStrength = 1;

Set.lambdaS1 = 5;
Set.lambdaS2 = 0.01*Set.lambdaS1;
Set.VTK = false;