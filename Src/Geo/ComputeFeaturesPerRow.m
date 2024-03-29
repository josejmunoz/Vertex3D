function row_features = ComputeFeaturesPerRow(Geo, cellsToAblate, row_features)
    %% Compute row analysis as if it were a wound
    nodesOfTheWound = getNodeNeighbours(Geo, cellsToAblate);
    pastCellsOfTheWound = cellsToAblate;
    cellsOfTheWound = horzcat(cellsToAblate, nodesOfTheWound(IsCell(Geo, nodesOfTheWound))');
    numRow = 1;
    while ~any(IsBorderCell(Geo, cellsOfTheWound))
        newWound_features = ComputeWoundFeatures(Geo, setdiff(cellsOfTheWound, pastCellsOfTheWound));
        for newField = fieldnames(newWound_features)'
            row_features.(strcat('Row', num2str(numRow), '_', newField{:})) = newWound_features.(newField{:});
        end

        pastCellsOfTheWound = cellsOfTheWound;
        nodesOfTheWound = getNodeNeighbours(Geo, cellsOfTheWound);
        cellsOfTheWound = horzcat(cellsToAblate, nodesOfTheWound(IsCell(Geo, nodesOfTheWound))');
        numRow = numRow + 1;
    end
end