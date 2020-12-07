function [g,K,Cell,Energy] = KgContractility(Cell,Y,Set)
%KGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here
%   g: is a vector
%   K: is a matrix

%% Set parameters
totalCells = Cell.n;

%% Initialize
dimensionsG = 1;%Set.NumTotalV*3;

K = zeros(dimensionsG,1); % Local cell residual

EnergyContractility = 0;

C = Set.cContractility;

Energy = 0;

[g] = computeGContractility(Cell, C);


end

function [kgContractility] = computeKGContractility(Cell, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

[Cell] = Cell.computeEdgeLengths();
Cell.EdgeLengths = Cell.EdgeLengths;

end

function [gContractility] = computeGContractility(Cell, C)
%COMPUTEGCONTRACTILITY Summary of this function goes here
%   Detailed explanation goes here

[Cell] = Cell.computeEdgeLengths();
Cell.EdgeLengths = Cell.EdgeLengths;

gContractility = 0.5 * (edgeLength - (edgeLength0 * C))^2;

end
