function [Geo] = CombineTwoGhostNodes(Geo, nodesToCombine)
%COMBINETWOGHOSTNODES Summary of this function goes here
%   Detailed explanation goes here

    if isempty([Geo.Cells(nodesToCombine).AliveStatus]) %% All of them need to be ghost nodes
        CellsToCombine = [Geo.Cells(nodesToCombine)];

        newCell = CellsToCombine(1);
        
        newCell.X = mean(vertcat(CellsToCombine.X));
        newCell.T = vertcat(CellsToCombine.T);
        
        %Replace old for new ID
        replacedTets = newCell.T(any(ismember(newCell.T, CellsToCombine(2).ID), 2), :);
        newCell.T(ismember(newCell.T, CellsToCombine(2).ID)) = newCell.ID;
        [~, nonRepeatIDs] = unique(sort(newCell.T, 2), 'rows');
        removedTets = newCell.T(setdiff(1:size(newCell.T, 1), nonRepeatIDs), :);
        newCell.T = newCell.T(nonRepeatIDs, :);
        
        for numCell = [Geo.Cells.ID]
            numCell
        end
    end
end

