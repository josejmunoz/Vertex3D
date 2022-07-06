function [Geo, Tnew, removedTets, replacedTets] = CombineTwoGhostNodes(Geo, Set, nodesToCombine)
%COMBINETWOGHOSTNODES Summary of this function goes here
%   Detailed explanation goes here
    Tnew = [];
    if isempty([Geo.Cells(nodesToCombine).AliveStatus]) %% All of them need to be ghost nodes
        CellsToCombine = [Geo.Cells(nodesToCombine)];

        newCell = CellsToCombine(1);
        
        newCell.X = mean(vertcat(CellsToCombine.X));
        newCell.T = vertcat(CellsToCombine.T);
        
        %Replace old for new ID
        replacingTets = ismember(newCell.T, CellsToCombine(2).ID);
        replacedTets = newCell.T(any(replacingTets, 2), :);
        newCell.T(replacingTets) = newCell.ID;
        %Remove repeated tets after replacement with new IDs
        [~, nonRepeatIDs] = unique(sort(newCell.T, 2), 'rows');
        removedTets = newCell.T(setdiff(1:size(newCell.T, 1), nonRepeatIDs), :);
        Tnew = newCell.T(any(replacingTets(sort(nonRepeatIDs), :), 2), :);
        newCell.T = newCell.T(nonRepeatIDs, :);
        % Removing Tets with the new cell twice or more within the Tet
        newCell.T(sum(ismember(newCell.T, nodesToCombine(1)), 2) > 1, :) = [];
        Tnew(sum(ismember(Tnew, nodesToCombine(1)), 2) > 1, :) = [];
        
        for numCell = [Geo.Cells.ID]
            currentTets = Geo.Cells(numCell).T;
            if any(any(ismember(currentTets, nodesToCombine(2))))
                %% Replace 'old' cell by 'new' cell
                replacingTets = ismember(currentTets, nodesToCombine(2));
                Geo.Cells(numCell).T(replacingTets) = nodesToCombine(1);
                currentTets = Geo.Cells(numCell).T;
                
                % Recalculate Ys
%                 if ~isempty(Geo.Cells(numCell).Y)
%                     for idTet = find(any(replacingTets, 2))'
%                         Geo.Cells(numCell).Y(idTet, :) = ComputeY(vertcat(Geo.Cells(currentTets(idTet, :)).X), newCell.X, length([Geo.Cells(currentTets(idTet, :)).AliveStatus]) > 1, Set);
%                     end
%                 end
                
                %% Remove repeated Tets
                % IDs are not ordered in the same way for different cells
                % but same Tet
                checkRepeatedTets = ismember(sort(currentTets, 2), sort(removedTets, 2), 'rows');
                checkRepatedCells = sum(ismember(currentTets, nodesToCombine(1)), 2) > 1;
                Geo.Cells(numCell).T(checkRepeatedTets | checkRepatedCells, :) = [];
                if ~isempty(Geo.Cells(numCell).Y)
                    Geo.Cells(numCell).Y(checkRepeatedTets | checkRepatedCells, :) = [];
                end
            end
        end
        
        %Update the 'new' cell
        Geo.Cells(nodesToCombine(1)) = newCell;
        %Remove the 'old' cell
        %Geo.Cells(nodesToCombine(2)).X = []; 
        Geo.Cells(nodesToCombine(2)).T = [];
    end
end

