function [Cell, Y] = simpleRemodelling(Cell, Y, tetrahedra, Set)
%SIMPLEREMODELLING Summary of this function goes here
%   Detailed explanation goes here

    %% Identify if any edge is shorter than it should
    remodellingCells = [];
    for numCell = find(Cell.DebrisCells == 0)'
        currentEdgeVertices = Cell.Cv{numCell};
        
        neighbours3D = unique(tetrahedra(any(ismember(tetrahedra, numCell), 2), :));
        neighbours3D = neighbours3D(ismember(neighbours3D, Cell.Int));
        for neighbourCell = neighbours3D'
            neighbourWoundEdges = vertcat(Cell.Cv{neighbourCell});
            idShareEdges = ismember(sort(currentEdgeVertices, 2), sort(neighbourWoundEdges, 2), 'rows');
            %% Apical edge
            apicalEdgeLength(end+1) = sum(Cell.EdgeLengths{numCell}(idShareEdges & Cell.EdgeLocation{numCell} == 3));
            %% Basal edge
            basalEdgeLength(end+1) = sum(Cell.EdgeLengths{numCell}(idShareEdges & Cell.EdgeLocation{numCell} == 2));
        end
    end
    
    if empty(remodellingCells) == 0
       % Intercalate whole cell in 3D: both cells do not share a face anymore 
        
    end
end

