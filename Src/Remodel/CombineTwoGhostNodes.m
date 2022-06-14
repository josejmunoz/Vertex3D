function [outputArg1,outputArg2] = CombineTwoGhostNodes(Geo, nodesToCombine)
%COMBINETWOGHOSTNODES Summary of this function goes here
%   Detailed explanation goes here

    if isempty([Geo.Cells(nodesToCombine).AliveStatus]) %% All of them need to be ghost nodes
        CellsToCombine = [Geo.Cells(nodesToCombine)];

        newCell = CellsToCombine(1);
    end
end

