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


end

