Geo.nx = 3; 
Geo.ny = 3; 
Geo.nz = 1;

Geo.cellsToAblate = [5];
 
Set.lambdaV = 10; 
 
Set.tend=300; 
Set.Nincr=300; 
 
Set.lambdaS1 = 0.2; 
Set.lambdaS2 = 0.2;  
 
Set.Substrate  = 0; 
Set.SubstrateZ = -0.5;

Set.Contractility = true;
Set.cLineTension = 0.001;

%Set.Remodelling = 0;
 
Set.OutputFolder='Result/Substrate';