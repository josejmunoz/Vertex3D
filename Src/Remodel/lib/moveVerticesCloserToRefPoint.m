function [Geo, Geo_n] = moveVerticesCloserToRefPoint(Geo, Geo_n, closeToNewPoint, cellNodesShared, cellToSplitFrom, ghostNode, Set)
%% Vertices connecting the two intercalating cells should be closer
    allT = vertcat(Geo.Cells.T);
    if ismember(ghostNode, Geo.XgBottom)
        allT_filtered = allT(any(ismember(allT, Geo.XgBottom), 2), :);
    elseif ismember(ghostNode, Geo.XgTop)
        allT_filtered = allT(any(ismember(allT, Geo.XgTop), 2), :);
    end

    %% Obtain point of reference of cells (i.e. 3 cell nodes, 1 ghost node)
    possibleRefTets = allT_filtered(sum(ismember(allT_filtered, cellNodesShared), 2) == 3, :);
    possibleRefTets = unique(sort(possibleRefTets(sum(ismember(possibleRefTets, cellNodesShared), 2) == 3, :), 2), 'rows');
    refTet = any(ismember(possibleRefTets, cellToSplitFrom), 2);
    refPoint_closer = Geo.Cells(cellToSplitFrom).Y(ismember(sort(Geo.Cells(cellToSplitFrom).T, 2), possibleRefTets(refTet, :), 'rows'), :);

    if sum(refTet) > 1
        error('moveVerticesCloserToRefPoint_line17');
    end

    %% Obtain vertices to change
    id_cellsToChange = setdiff(cellNodesShared, Geo.XgID);
    id_cellsToChange = id_cellsToChange([Geo.Cells(id_cellsToChange).AliveStatus] == 1);
    verticesToChange = [];
    for numCell = id_cellsToChange'
        news = sum(ismember(Geo.Cells(numCell).T, Geo.XgID), 2) > 2;
        news(sum(ismember(Geo.Cells(numCell).T, id_cellsToChange), 2) == 2 & sum(ismember(Geo.Cells(numCell).T, Geo.XgID), 2) == 2) = 1;
        news(sum(ismember(Geo.Cells(numCell).T, id_cellsToChange), 2) >= 3) = 1;
    
        % Remove only the tets from the domain it is not changing
        if ismember(ghostNode, Geo.XgBottom)
            news(any(ismember(Geo.Cells(numCell).T, Geo.XgTop), 2)) = 0;
        else
            news(any(ismember(Geo.Cells(numCell).T, Geo.XgBottom), 2)) = 0;
        end
        verticesToChange = vertcat(verticesToChange, Geo.Cells(numCell).T(news, :));
    end

    verticesToChange = unique(sort(verticesToChange, 2), 'rows');

    verticesToChange(~any(ismember(verticesToChange, Geo.XgID), 2), :) = [];

    % Remove debris cells (those won't move)
    idCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
    deadCells = idCells([Geo.Cells(idCells).AliveStatus] == 0);
    %verticesToChange(ismember(verticesToChange, deadCells), :) = [];

    %% cells that were splitted need to get closer
    % Cells that were joined need to get further
    cellsToGetFurther = intersect(possibleRefTets(1, :), possibleRefTets(2, :));
    cellsToGetCloser = setdiff(cellNodesShared, cellsToGetFurther);
    cellsToGetCloser = cellsToGetCloser([Geo.Cells(cellsToGetCloser).AliveStatus] == 1);

    refPoint_further = mean(Geo.Cells(cellsToGetFurther(1)).Y(ismember(sort(Geo.Cells(cellsToGetFurther(1)).T, 2), possibleRefTets, 'rows'), :), 1);
    farFromNewPoint = closeToNewPoint;

    [movementByCell{cellNodesShared}] = deal([]);

    for tetToCheck = verticesToChange'
        if any(ismember(tetToCheck, cellsToGetCloser))
            gettingCloser = 1;
        else
            gettingCloser = 0;
        end
        for nodeInTet = tetToCheck'
            if ~ismember(nodeInTet, Geo.XgID)
                newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
                
                if gettingCloser
                    avgPoint = refPoint_closer*(1-closeToNewPoint) + newPoint*closeToNewPoint;
                else
                    avgPoint = refPoint_further*(1-farFromNewPoint) + newPoint*farFromNewPoint;
                end

                if sum(ismember(tetToCheck, Geo.XgID)) == 3
                    movementByCell{nodeInTet} = vertcat(newPoint - avgPoint, movementByCell{nodeInTet});
                else
                    Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = refPoint_closer*(1-closeToNewPoint) + newPoint*closeToNewPoint;
                end
            end
        end
    end

    movementByCell_avg = cellfun(@(x) median(x, 1), movementByCell, 'UniformOutput', false);

    vertices_only1Cell = verticesToChange(sum(ismember(verticesToChange, Geo.XgID), 2) == 3, :);

    %% MOVING CELLS TO ACCOMODATE
    % Should be a gradient where closer cells move more and further cell
    % less
    for idCell = cellNodesShared'
        if any(ismember(idCell, cellsToGetCloser))
            gettingCloser = 1;
        else
            gettingCloser = 0;
        end
        currentVertices = vertices_only1Cell(any(ismember(vertices_only1Cell, idCell), 2), :);
        
        avgPoint = movementByCell_avg{idCell};

        if ~isempty(currentVertices)
            for tetToCheck = currentVertices'
                newPoint = Geo.Cells(idCell).Y(ismember(sort(Geo.Cells(idCell).T, 2), tetToCheck', 'rows'), :);
                if gettingCloser
                    Geo.Cells(idCell).Y(ismember(sort(Geo.Cells(idCell).T, 2), tetToCheck', 'rows'), :) = newPoint - avgPoint;
                else
                    Geo.Cells(idCell).Y(ismember(sort(Geo.Cells(idCell).T, 2), tetToCheck', 'rows'), :) = newPoint + avgPoint;
                end
            end
        end
    end

    %% Also the vertex middle Scutoid vertex
    for currentCell = cellNodesShared'
        middleVertexTet = all(ismember(Geo.Cells(currentCell).T, cellNodesShared), 2);
        Geo.Cells(currentCell).Y(middleVertexTet, :) = refPoint_closer*(1-closeToNewPoint) + Geo.Cells(currentCell).Y(middleVertexTet, :)*(closeToNewPoint);
    end

    %% Rebuild
    Geo = BuildXFromY(Geo_n, Geo);
    % Recalculating face centres here based on the previous
    % change
    Geo = Rebuild(Geo, Set);
    Geo = BuildGlobalIds(Geo);
    Geo   = UpdateMeasures(Geo);
    Geo_n = Geo;
end