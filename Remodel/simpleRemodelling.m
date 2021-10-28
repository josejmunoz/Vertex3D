function [Cell, Y, tetrahedra] = simpleRemodelling(Cell, Y0, Yn, Y, CellInput, tetrahedra_, Tetrahedra_weights, X, X_IDs, SCn, XgID, Cn, Set)
%SIMPLEREMODELLING Summary of this function goes here
%   Detailed explanation goes here

    tetrahedra = tetrahedra_.DataRow;

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
            apicalEdgeLength(end+1) = apicalArea / sum(Cell.EdgeLengths{numCell}(idShareEdges & Cell.EdgeLocation{numCell} == 3));
            %% Basal edge
            basalEdgeLength(end+1) = basalArea / sum(Cell.EdgeLengths{numCell}(idShareEdges & Cell.EdgeLocation{numCell} == 2));
            
            if basalEdgeLength(end) > Set.MinEdgeLength || apicalEdgeLength(end) > Set.MinEdgeLength
                remodellingCells(end+1, 1:2) = [numCell, neighbourCell];
            end
        end
        apicalLengths{numCell} = apicalEdgeLength;
        basalLengths{numCell} = basalEdgeLength;
    end
    
    if isempty(remodellingCells) == 0
       %% Intercalate whole cell in 3D: both cells do not share a face anymore 
       for numIntercalation = 1:size(remodellingCells, 1)
           currentIntercalation = remodellingCells(numIntercalation, :); % Face to remove
           tetsToChange = tetrahedra(sum(ismember(tetrahedra, currentIntercalation), 2) >=2, :);
           verticesToChange = unique(tetsToChange);
           nodesToChange = verticesToChange(ismember(verticesToChange, Cell.Int));
           tetsToChangeAll_IDs = sum(ismember(tetrahedra, verticesToChange), 2) >= 3 ;
           tetsToChangeAll = tetrahedra(tetsToChangeAll_IDs , :);
           
           newConnectedNodes = nodesToChange(ismember(nodesToChange, currentIntercalation) == 0);
           
           %% Reconstruct Tets
           % New triangles
           trianglesConnectivity = nchoosek(nodesToChange, 3);
           trianglesConnectivity(sum(ismember(trianglesConnectivity, currentIntercalation), 2) >=2, :) = [];
           
           % Relationships: 1 ghost node, three cell nodes
           tetsToChange_1 = tetrahedra(sum(ismember(tetrahedra, nodesToChange), 2) >= 3 , :);
           Twg_vertices_1 = horzcat(repmat(trianglesConnectivity, 2, 1), tetsToChange_1(:, 4));
           
           % Relationships: 2 ghost nodes, two cell nodes
           tetsToChange_2 = tetsToChangeAll(sum(ismember(tetsToChangeAll, nodesToChange), 2) == 2 , :);
           Twg_vertices_2 = tetsToChange_2(sum(ismember(tetsToChange_2, currentIntercalation), 2) == 2, :);
           Twg_vertices_2(Twg_vertices_2 == currentIntercalation(1)) = newConnectedNodes(1);
           Twg_vertices_2(Twg_vertices_2 == currentIntercalation(2)) = newConnectedNodes(2);
           tetsToChange_2(sum(ismember(tetsToChange_2, currentIntercalation), 2) == 2, :) = Twg_vertices_2;
           Twg_vertices_2 = tetsToChange_2;
           
           % Relationships: 1 cell node and 3 ghost nodes
           Twg_vertices_3_1 = horzcat(newConnectedNodes, X_IDs.topFaceIds(newConnectedNodes)', repmat(tetsToChange_1(1:2, 4)', 2, 1));
           Twg_vertices_3_2 = horzcat(newConnectedNodes, X_IDs.bottomFaceIds(newConnectedNodes)', repmat(tetsToChange_1(3:4, 4)', 2, 1));
           
           %% New tetrahedra substitution
           tetrahedra(tetsToChangeAll_IDs, :) = vertcat(Twg_vertices_1(1:2, :), Twg_vertices_2(1:size(Twg_vertices_2, 1)/2, :), Twg_vertices_3_1, ...
              Twg_vertices_1(3:4, :), Twg_vertices_2(size(Twg_vertices_2, 1)/2+1:end, :), Twg_vertices_3_2);
           
           %% New Y_s
           for numTetrahedron = find(tetsToChangeAll_IDs)'
               Y.DataRow(numTetrahedron, :) = mean(X(tetrahedra(numTetrahedron, :), :)) ./ Tetrahedra_weights(numTetrahedron, :);
           end
           
           %% Remove face
           faceToRemove = find(all(ismember(Cell.AllFaces.Nodes, currentIntercalation), 2));
           Cell.AllFaces=Cell.AllFaces.Remove(faceToRemove);
           SCn=SCn.Remove(faceToRemove);
           Cell.FaceCentres=Cell.FaceCentres.Remove(faceToRemove);
           
           % Remove face on cells
           for numCell = nodesToChange'
               idsToRemove = ismember(Cell.Faces{numCell}.FaceCentresID, faceToRemove);
               Cell.Faces{numCell}.FaceCentresID(idsToRemove) = [];
               Cell.Faces{numCell}.Vertices(idsToRemove) = [];
               Cell.Faces{numCell}.Tris(idsToRemove) = [];
               Cell.Faces{numCell}.nFaces = Cell.Faces{numCell}.nFaces - sum(idsToRemove);
               Cell.cNodes{numCell}(idsToRemove) = [];
           end
           
           %% Add new face
           Cell.cNodes{nodesToChange(1)}(end) = nodesToChange(2);
           Cell.cNodes{nodesToChange(2)}(end) = nodesToChange(1);
           
           %% Rebuild cells
           Cell.AssembleNodes=nodesToChange;
           [Cell, nC, SCn, flag]=ReBuildCells(Cell, tetrahedra_, Y, X, SCn);
           
           if ~flag
               Cell.AllFaces=Cell.AllFaces.CheckInteriorFaces(XgID);
               [Cell]=ComputeCellVolume(Cell,Y);
               Cell = Cell.computeEdgeLengths(Y);
               for jj=1:Cell.n
                   Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
                   Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
               end
           end
           
           %% Solve modelling step with only those vertices
           [Dofs] = GetDOFs(Y, Cell, Set, isempty(Set.InputSegmentedImage) == 0);
           [Dofs] = updateRemodelingDOFs(Dofs, find(tetsToChangeAll_IDs), nC, Y);
           
           Cell.RemodelledVertices=[find(tetsToChangeAll_IDs)', nC+Y.n];
           if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end

           [Cell,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);
       end
    end
end

