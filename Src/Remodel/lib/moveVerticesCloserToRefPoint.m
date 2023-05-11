function [Geo, Geo_n] = moveVerticesCloserToRefPoint(Geo, Geo_n, cellNodesShared, cellToSplitFrom, Set)
%% Vertices connecting the two intercalating cells should be closer
    
    closeToNewPoint = 0.1;

    allT = vertcat(Geo.Cells.T);
    if ismember(ghostNode, Geo.XgBottom)
        allT_filtered = allT(any(ismember(allT, Geo.XgBottom), 2), :);
    elseif ismember(ghostNode, Geo.XgTop)
        allT_filtered = allT(any(ismember(allT, Geo.XgTop), 2), :);
    end

    % Vertices of cells (i.e. 3 cell nodes, 1 ghost node)
    verticesToChange = allT_filtered(sum(ismember(allT_filtered, cellNodesShared), 2) == 3, :);
    verticesToChange = unique(sort(verticesToChange(sum(ismember(verticesToChange, cellNodesShared), 2) == 3, :), 2), 'rows');

    refTet = any(ismember(verticesToChange, cellToSplitFrom), 2);
    refPoint = Geo.Cells(cellToSplitFrom).Y(ismember(sort(Geo.Cells(cellToSplitFrom).T, 2), verticesToChange(refTet, :), 'rows'), :);

    if sum(refTet) > 1
        disp('error');
    end

    cellsConnected = intersect(verticesToChange(1, :), verticesToChange(2, :));

    verticesToChange(refTet, :) = [];

    middleVertexToChange = allT_filtered(sum(ismember(allT_filtered, cellsConnected), 2) == 2 & sum(ismember(allT_filtered, Geo.XgID), 2) == 2, :);
    middleVertexToChange = unique(sort(middleVertexToChange, 2), 'rows');

    verticesToChange = vertcat(verticesToChange, middleVertexToChange);

    for tetToCheck = verticesToChange'
        for nodeInTet = tetToCheck'
            if ~ismember(nodeInTet, Geo.XgID)
                newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);

                Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = refPoint*(1-closeToNewPoint) + newPoint*closeToNewPoint;
            end
        end
    end

    %closeToNewPoint = 0.5;
    % Also the vertex middle Scutoid vertex
    for currentCell = cellNodesShared'
        middleVertexTet = all(ismember(Geo.Cells(currentCell).T, cellNodesShared), 2);
        Geo.Cells(currentCell).Y(middleVertexTet, :) = refPoint*(1-closeToNewPoint) + Geo.Cells(currentCell).Y(middleVertexTet, :)*(closeToNewPoint);
    end

    Geo = BuildXFromY(Geo_n, Geo);
    % Recalculating face centres here based on the previous
    % change
    Geo = Rebuild(Geo, Set);
    Geo = BuildGlobalIds(Geo);
    Geo   = UpdateMeasures(Geo);
    Geo_n = Geo;
end