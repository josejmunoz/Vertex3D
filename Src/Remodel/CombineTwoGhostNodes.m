function [Geo] = CombineTwoGhostNodes(Geo, nodesToCombine)
%COMBINETWOGHOSTNODES Summary of this function goes here
%   Detailed explanation goes here

    if isempty([Geo.Cells(nodesToCombine).AliveStatus]) %% All of them need to be ghost nodes
        CellsToCombine = [Geo.Cells(nodesToCombine)];

        newCell = CellsToCombine(1);
        
        newCell.X = mean(vertcat(CellsToCombine.X));
        newCell.T = vertcat(CellsToCombine.T);
        
        %Replace old for new ID
        newCell.T(ismember(newCell.T, CellsToCombine(2).ID)) = newCell.ID;
        [~, nonRepeatIDs] = unique(sort(newCell.T, 2), 'rows');
        removedTets = newCell.T(setdiff(1:size(newCell.T, 1), nonRepeatIDs), :);
        newCell.T = newCell.T(nonRepeatIDs, :);
        
        for numCell = [Geo.Cells.ID]
            currentTets = Geo.Cells(numCell).T;
            if any(any(ismember(currentTets, nodesToCombine(2))))
                %% Replace 'old' cell by 'new' cell
                replacedTets = ismember(currentTets, nodesToCombine(2));
                Geo.Cells(numCell).T(replacedTets) = nodesToCombine(1);
                
                % Recalculate Ys
                for idTet = find(any(replacedTets, 2))
                    Geo.Cells(numCell).Y(idTet, :) = ComputeY(vertcat(Geo.Cells(currentTets(idTet, :).X), newCell.X, -1, false);
                end
                
                % IDs are not ordered in the same way for different cells
                % but same Tet
                checkRepeatedTets = ismember(sort(currentTets, 2), sort(removedTets, 2), 'rows');
                Geo.Cells(numCell).T(checkRepeatedTets, :) = [];
                Geo.Cells(numCell).Y(checkRepeatedTets, :) = [];
                

            end
        end
        
        %Update the 'new' cell
        Geo.Cells(nodesToCombine(1)) = newCell;
        %Remove the 'old' cell
        Geo.Cells(nodesToCombine(2)).X = []; 
        Geo.Cells(nodesToCombine(2)).T = []; 
    end
end

