function [] = gContractility(Cell,Y,Set)
%GCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

%% Set parameters
currentCell=Cell.n;

%% Initialize
dimensionsG = %Set.NumTotalV*3;

g = zeros(dimensionsG,1); % Local cell residual

EnergyContractility=0;

C=Set.cContractility;

end


function [gContractility] = computeGContractility(C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here



end