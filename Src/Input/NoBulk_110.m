Set.InputGeo = 'VertexModelTime';
% 40 cells; 3 cells to ablate
% 110 cells; 7 cells to ablate
Set.TotalCells = 110; % 197 cell crashes
Geo.cellsToAblate = 1:3;

Set.InPlaneElasticity = false;
Set.mu_bulk	= 0; %% 1000
Set.lambda_bulk	= 0; %% 500

Set.nu	= 500;
Set.Nincr = 61 * 100;

Set.lambdaB	= 5; % Smaller number, more energy (resulting in bigger triangles)
Set.lambdaR	= 0.3; % Aspect ratio of triangles 

Set.lambdaV = 1;

Set.kSubstrate = 1;

Set.cLineTension = 0.01;
Set.purseStringStrength = 0; % 11 for an hour looks good % 25 for tests
Set.noiseContractility = 0;
Set.DelayedAdditionalContractility = 0;

% Soft < 0
% Stiff > 0
Set.Remodelling = false;
Set.RemodelStiffness = 1;

Set.lambdaS1 = 10;
Set.lambdaS2 = Set.lambdaS1/10;
Set.lambdaS3 = Set.lambdaS1/10;
Set.lambdaS4 = Set.lambdaS3;

Set.VTK = true;