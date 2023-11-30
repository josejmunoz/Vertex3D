Set.InputGeo = 'Bubbles_Cyst';
Set.TotalCells = 20;
Set.s = 1.5*10;
Set.f = 0.5*10;

Set.Ablation = false;

Set.InPlaneElasticity = false;
Set.nu	= 1;
Set.Nincr = 61;

Set.lambdaB	= 1; % Smaller number, more energy
Set.lambdaR	= 0.3; 

Set.lambdaV = 1000;

Set.Substrate = 0;
Set.kSubstrate = 0;

Set.cLineTension = 0.1;

Set.lambdaS1 = 8;
Set.lambdaS2 = 0.01*Set.lambdaS1;
Set.lambdaS3 = Set.lambdaS2;
Set.lambdaS4 = Set.lambdaS2;

Set.Remodelling = 1;
Set.RemodelStiffness = 0.5;
Set.VTK = true;