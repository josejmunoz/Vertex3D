function [Geo, Tnew, Ynew, removedTets, replacedTets] = CombineTwoNodes(Geo, Set, nodesToCombine)
%COMBINETWOGHOSTNODES Summary of this function goes here
%   Detailed explanation goes here

    oldGeo = Geo;
    CellsToCombine = [Geo.Cells(nodesToCombine)];
    
    nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    newCell = CellsToCombine(1);
    
    oldYs = Geo.Cells(nodesToCombine(2)).Y;
    oldTs = Geo.Cells(nodesToCombine(2)).T;
    
    newCell.X = mean(vertcat(CellsToCombine.X));
    newCell.T = vertcat(CellsToCombine.T);
    newCell.Y = vertcat(CellsToCombine.Y);
    
    %Replace old for new ID
    replacingTets = ismember(newCell.T, CellsToCombine(2).ID);
    replacedTets = newCell.T(any(replacingTets, 2), :);
    newCell.T(replacingTets) = newCell.ID;
    %Remove repeated tets after replacement with new IDs
    [~, nonRepeatIDs] = unique(sort(newCell.T, 2), 'rows');
    removedTets = newCell.T(setdiff(1:size(newCell.T, 1), nonRepeatIDs), :);
    newCell.Y = newCell.Y(nonRepeatIDs, :);
    newCell.T = newCell.T(nonRepeatIDs, :);
    % Removing Tets with the new cell twice or more within the Tet
    newCell.Y(sum(ismember(newCell.T, nodesToCombine(1)), 2) > 1, :) = [];
    newCell.T(sum(ismember(newCell.T, nodesToCombine(1)), 2) > 1, :) = [];
    
    for numCell = [Geo.Cells.ID]
        currentTets = Geo.Cells(numCell).T;
        if any(any(ismember(currentTets, [nodesToCombine(2) nodesToCombine(1)])))
            %% Replace 'old' cell by 'new' cell
            replacingTets = ismember(currentTets, nodesToCombine(2));
            Geo.Cells(numCell).T(replacingTets) = nodesToCombine(1);
            currentTets = Geo.Cells(numCell).T;
            
            %% Remove repeated Tets
            % IDs are not ordered in the same way for different cells
            % but same Tet
            checkRepeatedTets = ismember(sort(currentTets, 2), sort(removedTets, 2), 'rows');
            checkRepatedCells = sum(ismember(currentTets, nodesToCombine(1)), 2) > 1;
            Geo.Cells(numCell).T(checkRepeatedTets | checkRepatedCells, :) = [];
            if ~ismember(numCell, Geo.XgID)
                Geo.Cells(numCell).Y(checkRepeatedTets | checkRepatedCells, :) = [];
            end
            replacingTets(checkRepeatedTets | checkRepatedCells, :) = [];
        end
    end
    
    %Update the 'new' cell
    Geo.Cells(nodesToCombine(1)) = newCell;
    %Remove the 'old' cell
    Geo.Cells(nodesToCombine(2)).T = [];
    Geo.Cells(nodesToCombine(2)).Y = [];
    Geo.Cells(nodesToCombine(2)).AliveStatus = [];
    Geo.Cells(nodesToCombine(2)).Area = [];
    Geo.Cells(nodesToCombine(2)).Area0 = [];
    Geo.Cells(nodesToCombine(2)).Vol = [];
    Geo.Cells(nodesToCombine(2)).Vol0 = [];
    Geo.Cells(nodesToCombine(2)).Faces = [];
    Geo.Cells(nodesToCombine(2)).cglobalIds = [];
    Geo.Cells(nodesToCombine(2)).globalIds = [];
    Geo.Cells(nodesToCombine(2)).ExternalLambda = [];
    Geo.Cells(nodesToCombine(2)).InternalLambda = [];
    Geo.Cells(nodesToCombine(2)).SubstrateLambda = [];
    
    Tnew = newCell.T;
    Ynew = newCell.Y;

end

