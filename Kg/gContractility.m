function [EnergyContractility] = gContractility(Cell,Y,Set)
%GCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

%% Set parameters
totalCells = Cell.n;

%% Initialize
dimensionsG = 1;%Set.NumTotalV*3;

g = zeros(dimensionsG,1); % Local cell residual

EnergyContractility = 0;

C = Set.cContractility;

[gContractility] = computeGContractility(Cell, C);

end


function [gContractility] = computeGContractility(Cell, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

[Cell] = Cell.computeEdgeLengths();
Cell.EdgeLengths = Cell.EdgeLengths;

gContractility = 0.5 * (edgeLength - (edgeLength0 * C))^2;

end