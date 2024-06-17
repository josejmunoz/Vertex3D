Set.InputGeo = 'VertexModelTime';
% 40 cells; 3 cells to ablate
% 110 cells; 7 cells to ablate
Set.TotalCells = 110; % 197 cell crashes
Geo.cellsToAblate = 1:7;

Set.InPlaneElasticity = false;
Set.mu_bulk	= 0; %% 1000
Set.lambda_bulk	= 0; %% 500

Set.nu	= 5000;
Set.Nincr = 61;

Set.lambdaB	= 1; % Smaller number, more energy (resulting in bigger triangles)
Set.lambdaR	= 0.3; % Aspect ratio of triangles 

Set.lambdaV = 10;

Set.kSubstrate = 100;

Set.cLineTension = 0.05;
Set.purseStringStrength = 12; % 11 for an hour looks good % 25 for tests
Set.noiseContractility = 0;
Set.DelayedAdditionalContractility = 0;

% Soft < 0
% Stiff > 0
Set.RemodelStiffness = 0;

Set.lambdaS1 = 10;
Set.lambdaS2 = Set.lambdaS1/10;
Set.lambdaS3 = Set.lambdaS1/10;
Set.lambdaS4 = Set.lambdaS3;

Set.VTK = true;