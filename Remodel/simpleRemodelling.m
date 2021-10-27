function [Cell, Y] = simpleRemodelling(Cell, Y, tetrahedra, Set)
%SIMPLEREMODELLING Summary of this function goes here
%   Detailed explanation goes here

    %% Identify if any edge is shorter than it should
    remodellingCells = [];
    apicalLengths = {};
    basalLengths = {};
    for numCell = find(Cell.DebrisCells == 0)'
        currentEdgeVertices = Cell.Cv{numCell};
        
        trianglesArea = Cell.SAreaTri{numCell};
        apicalFaceCentres = abs(Cell.ApicalVertices{numCell}(Cell.ApicalVertices{numCell} < 0));
        apicalTriangles = all(ismember(Cell.Tris{numCell}(:, 1:2), Cell.ApicalVertices{numCell}), 2) & ismember(Cell.Tris{numCell}(:, 3), apicalFaceCentres);
        apicalArea = sum(trianglesArea(apicalTriangles));
        
        basalFaceCentres = abs(Cell.BasalVertices{numCell}(Cell.BasalVertices{numCell} < 0));
        basalTriangles = all(ismember(Cell.Tris{numCell}(:, 1:2), Cell.BasalVertices{numCell}), 2) & ismember(Cell.Tris{numCell}(:, 3), basalFaceCentres);
        basalArea = sum(trianglesArea(basalTriangles));
        
        neighbours3D = unique(tetrahedra(any(ismember(tetrahedra, numCell), 2), :));
        neighbours3D = neighbours3D(ismember(neighbours3D, Cell.Int));
        apicalEdgeLength = [];
        basalEdgeLength = [];
        for neighbourCell = neighbours3D'
            if neighbourCell == numCell
                continue
            end
            neighbourWoundEdges = vertcat(Cell.Cv{neighbourCell});
            idShareEdges = ismember(sort(currentEdgeVertices, 2), sort(neighbourWoundEdges, 2), 'rows');
            %% Apical edge
            apicalEdgeLength(end+1) = sum(Cell.EdgeLengths{numCell}(idShareEdges & Cell.EdgeLocation{numCell} == 3)) / apicalArea;
            %% Basal edge
            basalEdgeLength(end+1) = sum(Cell.EdgeLengths{numCell}(idShareEdges & Cell.EdgeLocation{numCell} == 2))  / basalArea;
            
            if basalEdgeLength(end) < Set.MinEdgeLength || apicalEdgeLength(end) < Set.MinEdgeLength
                remodellingCells(end+1, 1:2) = [numCell, neighbourCell];
            end
        end
        apicalLengths{numCell} = apicalEdgeLength;
        basalLengths{numCell} = basalEdgeLength;
    end
    
    if empty(remodellingCells) == 0
       % Intercalate whole cell in 3D: both cells do not share a face anymore 
        
    end
end

