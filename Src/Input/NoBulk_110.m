Set.InputGeo = 'VertexModelTime';
% 40 cells; 3 cells to ablate
% 110 cells; 7 cells to ablate
Set.TotalCells = 110;
Geo.cellsToAblate = 1:7;

Set.InPlaneElasticity = false;
Set.mu_bulk	= 0; %% 1000
Set.lambda_bulk	= 0; %% 500

Set.nu	= 5000;
Set.Nincr = 61*2;

Set.lambdaB	= 1; % Smaller number, more energy
Set.lambdaR	= 0.3; 

Set.lambdaV = 10;

Set.kSubstrate = 100;

Set.cLineTension = 0.05;
Set.purseStringStrength = 30;

Set.lambdaS1 = 20;
Set.lambdaS3 = Set.lambdaS1/10;
Set.VTK = true;