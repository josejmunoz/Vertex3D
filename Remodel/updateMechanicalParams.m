function [Set] = updateMechanicalParams(Set, initValues, maxSteps, numStep)
%UPDATEMECHANICALPARAMS Summary of this function goes here
%   Detailed explanation goes here
Set.lambdaV = initValues.lambdaV * (numStep / maxSteps);
Set.lambdaV_Debris = initValues.lambdaV_Debris * (numStep / maxSteps);
Set.lambdaS1 = initValues.lambdaS1 * (numStep / maxSteps);
Set.lambdaS2 = initValues.lambdaS2 * (numStep / maxSteps);
Set.lambdaS3 = initValues.lambdaS3 * (numStep / maxSteps);
Set.cLineTension = initValues.cLineTension * (numStep / maxSteps);
Set.mu_bulk = initValues.mu_bulk * (numStep / maxSteps);
Set.lambda_bulk = initValues.lambda_bulk * (numStep / maxSteps);
Set.kSubstrate = initValues.kSubstrate * (numStep / maxSteps);
Set.cPurseString = initValues.cPurseString * (numStep / maxSteps);
Set.cLateralCables = initValues.cLateralCables * (numStep / maxSteps);

Set.nu = initValues.nu * (maxSteps^2);
Set.nu0 = Set.nu;
end

