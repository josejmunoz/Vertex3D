function [remodellingCells, apicalLengths,basalLengths] = identifyEdgesToIntercalate(Cell, tetrahedra, Set)
%IDENTIFYEDGESTOINTERCALATE Summary of this function goes here
%   Detailed explanation goes here
    remodellingCells = [];
    apicalLengths = {};
    basalLengths = {};
    for numCell = find(Cell.DebrisCells == 0)'
        currentEdgeVertices = Cell.Cv{numCell};
        edgeLocation = Cell.EdgeLocation{numCell};
        edgeLengths = Cell.EdgeLengths{numCell};
        
        apicalPerimeterCell = sum(edgeLengths(edgeLocation == 3));
        basalPerimeterCell = sum(edgeLengths(edgeLocation == 2));
        
%         trianglesArea = Cell.SAreaTri{numCell};
%         apicalFaceCentres = abs(Cell.ApicalVertices{numCell}(Cell.ApicalVertices{numCell} < 0));
%         apicalTriangles = all(ismember(Cell.Tris{numCell}(:, 1:2), Cell.ApicalVertices{numCell}), 2) & ismember(Cell.Tris{numCell}(:, 3), apicalFaceCentres);
%         apicalArea = sum(trianglesArea(apicalTriangles));
        
%         basalFaceCentres = abs(Cell.BasalVertices{numCell}(Cell.BasalVertices{numCell} < 0));
%         basalTriangles = all(ismember(Cell.Tris{numCell}(:, 1:2), Cell.BasalVertices{numCell}), 2) & ismember(Cell.Tris{numCell}(:, 3), basalFaceCentres);
%         basalArea = sum(trianglesArea(basalTriangles));
        
        neighbours3D = unique(tetrahedra(any(ismember(tetrahedra, numCell), 2), :));
        neighbours3D = neighbours3D(ismember(neighbours3D, Cell.Int));
        neighbours3D(neighbours3D == numCell) = [];
        apicalEdgeLength = [];
        basalEdgeLength = [];
        for neighbourCell = neighbours3D'
            if neighbourCell == numCell
                continue
            end
            neighbourWoundEdges = vertcat(Cell.Cv{neighbourCell});
            idShareEdges = ismember(sort(currentEdgeVertices, 2), sort(neighbourWoundEdges, 2), 'rows');
            %% Apical edge
            apicalEdges = idShareEdges & Cell.EdgeLocation{numCell} == 3;
            adjacentApicalEdges = any(ismember(currentEdgeVertices, unique(currentEdgeVertices(apicalEdges, :))), 2) & Cell.EdgeLocation{numCell} == 3 & ~apicalEdges;
            adjacentApicalEdgesLength = sum(Cell.EdgeLengths{numCell}(adjacentApicalEdges));
            apicalEdgeLength(end+1) = apicalPerimeterCell / (sum(Cell.EdgeLengths{numCell}(apicalEdges)) / adjacentApicalEdgesLength^2);
            %% Basal edge
            basalEdges = idShareEdges & Cell.EdgeLocation{numCell} == 2;
            adjacentBasalEdges = any(ismember(currentEdgeVertices, unique(currentEdgeVertices(basalEdges, :))), 2) & Cell.EdgeLocation{numCell} == 2 & ~basalEdges;
            adjacentBasalEdgesLength = sum(Cell.EdgeLengths{numCell}(adjacentBasalEdges));
            basalEdgeLength(end+1) = basalPerimeterCell / (sum(Cell.EdgeLengths{numCell}(basalEdges)) / adjacentBasalEdgesLength^2);
            
            %basalEdgeLength(end) = basalEdgeLength(end) * 
            
            if (basalEdgeLength(end) > Set.MinEdgeLength || apicalEdgeLength(end) > Set.MinEdgeLength) 
                tetsToChange = tetrahedra(sum(ismember(tetrahedra, [numCell neighbourCell]), 2) >=2, :);
                verticesToChange = unique(tetsToChange);
                nodesToChange = verticesToChange(ismember(verticesToChange, Cell.Int));
                if sum(ismember(nodesToChange, find(Cell.DebrisCells))) <= 1
                    remodellingCells(end+1, 1:4) = [numCell, neighbourCell, apicalEdgeLength(end), basalEdgeLength(end)];
                end
            end
        end
        apicalLengths{numCell} = apicalEdgeLength;
        basalLengths{numCell} = basalEdgeLength;
    end
end

