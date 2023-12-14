Set.InputGeo = 'Bubbles_Cyst';
Set.typeOfEllipsoid = 'ellipsoid';
Set.ellipsoid_axis1 = 33/2;
Set.ellipsoid_axis2 = 31/2;
Set.ellipsoid_axis3 = 23/2;
Set.lumen_axis1 = 23/2;
Set.lumen_axis2 = 20/2;
Set.lumen_axis3 = 12/2;
Set.cell_V0 = 440;
Set.lumen_V0 = 4308;
Set.TotalCells = 30;
Set.s = 1.5*10;
Set.f = 0.5*10;

Set.Ablation = false;

Set.InPlaneElasticity = false;
Set.nu	= 1;
Set.Nincr = 61;

Set.lambdaB	= 1; % Smaller number, more energy
Set.lambdaR	= 0.3; 

Set.lambdaV = 100;

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