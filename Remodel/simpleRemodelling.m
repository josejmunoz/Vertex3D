function [Cell,Y,Yn,SCn,tetrahedra_,X,Dofs,Cn,Set] = simpleRemodelling(Cell, Y0, Yn, Y, CellInput, tetrahedra_, X, X_IDs, SCn, XgID, Cn, verticesInfo, neighboursNetwork, Dofs, Set)
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
           Yp = Y; tetrahedra_p = tetrahedra;
           currentIntercalation = remodellingCells(numIntercalation, :); % Face to remove
           tetsToChange = tetrahedra(sum(ismember(tetrahedra, currentIntercalation), 2) >=2, :);
           verticesToChange = unique(tetsToChange);
           nodesToChange = verticesToChange(ismember(verticesToChange, Cell.Int));
           tetsToChangeAll_IDs = sum(ismember(tetrahedra, verticesToChange), 2) >= 3 ;
           
           newConnectedNodes = nodesToChange(ismember(nodesToChange, currentIntercalation) == 0);
           
           % Reconstruct Tets
%            index = find(tetsToChangeAll_IDs)';
%            figure, tetramesh(tetrahedra(index, :), X);
           %% Option 1: New triangles
           % Change old connection for new connections
           trianglesConnectivity_all = verticesInfo.connectedCells;
           trianglesConnectivity_all(all(ismember(trianglesConnectivity_all, nodesToChange), 2), :) = sort([repmat(newConnectedNodes', 2, 1), currentIntercalation'], 2);
           neighboursNetwork(all(ismember(neighboursNetwork, currentIntercalation), 2), :) = [newConnectedNodes(1) newConnectedNodes(2); newConnectedNodes(2) newConnectedNodes(1)];
           
           % Remove vertices from edges
           vertices2DToChange = find(all(ismember(trianglesConnectivity_all, nodesToChange), 2));
           verticesConnectedToCells = [currentIntercalation(2) currentIntercalation(1)]; %% IT SHOULD BE DIFFERENT ORDER AS IN THE TRIANGULATION
           for numCellOriginal = 1:size(verticesInfo.edges, 1)
               currentEdges = verticesInfo.edges{numCellOriginal};
               if all(ismember(vertices2DToChange, currentEdges))
                   vertexToRemove = vertices2DToChange(ismember(verticesConnectedToCells, numCellOriginal));
                   edgesToRemoveIDs = find(any(currentEdges == vertexToRemove, 2));
                   edgesToRemove = currentEdges(edgesToRemoveIDs, :);
                   
                   if edgesToRemoveIDs(1) == 1 && edgesToRemoveIDs(2) == size(currentEdges, 1)
                       currentEdges = [edgesToRemove(2), edgesToRemove(3) ; currentEdges(edgesToRemoveIDs(1)+1:end-1, :)];
                   else
                       currentEdges = [currentEdges(1:edgesToRemoveIDs(1)-1, :); ...
                            edgesToRemove(1), edgesToRemove(end) ; currentEdges(edgesToRemoveIDs(2)+1:end, :)];
                   end
               end
               verticesInfo.edges{numCellOriginal} = currentEdges;
           end
           
           % Add vertices in edges corresponding to the new formed ones at
           % the newConnectedNodes
           for numCellOriginal = newConnectedNodes'
               currentEdges = verticesInfo.edges{numCellOriginal}(any(ismember(verticesInfo.edges{numCellOriginal}, vertices2DToChange), 2), :);
               vertexToAdd = vertices2DToChange(ismember(vertices2DToChange, currentEdges) == 0);
               vertexAdded = vertices2DToChange(ismember(vertices2DToChange, currentEdges));
               adjacentVertices = currentEdges(ismember(currentEdges, vertices2DToChange) == 0);
               
               adjacentVertexToVertexToAdd = adjacentVertices(sum(ismember(trianglesConnectivity_all(adjacentVertices, :), trianglesConnectivity_all(vertexToAdd, :)), 2)>1);
               edgeToSplit = find(all(ismember(verticesInfo.edges{numCellOriginal}, [adjacentVertexToVertexToAdd vertexAdded]), 2));
               
               % Select the adjacent vertex that shares 2 neighbours
               if adjacentVertexToVertexToAdd == adjacentVertices(1)
                   newEdges = [adjacentVertexToVertexToAdd vertexToAdd; vertexToAdd vertexAdded];
               else
                   newEdges = [vertexAdded vertexToAdd; vertexToAdd adjacentVertexToVertexToAdd]; 
               end
               verticesInfo.edges{numCellOriginal} = [verticesInfo.edges{numCellOriginal}(1:edgeToSplit-1, :); ...
                      newEdges; verticesInfo.edges{numCellOriginal}(edgeToSplit+1:end, :)];
           end
           verticesInfo.connectedCells = trianglesConnectivity_all;
           
           [Twg_bottom] = createTetrahedra(trianglesConnectivity_all, neighboursNetwork, verticesInfo.edges, Cell.Int', X_IDs.bottomFaceIds, X_IDs.bottomVerticesIds);
           [Twg_top] = createTetrahedra(trianglesConnectivity_all, neighboursNetwork, verticesInfo.edges, Cell.Int', X_IDs.topFaceIds, X_IDs.topVerticesIds);
           newTets = vertcat(Twg_top, Twg_bottom);
           newTets(all(ismember(newTets,XgID),2),:) = [];
           [changedTets, Locb] = ismember(tetrahedra, newTets, 'rows');
           [changedTets_2, Locb] = ismember(newTets, tetrahedra, 'rows');
           tetIds = 1:tetrahedra_.n;
           missingTets = tetIds(ismember(tetIds, Locb) == 0);
           newTetsModified = newTets(changedTets_2==0, :);
           
%            %%%%%% CHANGE ORDER OF NEW TETS TO MATCH THE OLD ONES
           tetrahedra(missingTets, :) = repmat([0 0 0 0], length(missingTets), 1);
           tetrahedra_.DataRow(missingTets, :) = repmat([0 0 0 0], length(missingTets), 1);
           Y.DataRow(missingTets, :) = repmat([-100 -100 -100], length(missingTets), 1);
           tetrahedra(tetrahedra_.n+1:tetrahedra_.n+size(newTetsModified, 1), :) = newTetsModified;

           % New nodes
           tetsToChange_1 = tetrahedra(sum(ismember(tetrahedra, nodesToChange), 2) >= 3 , :);
           for numX = 1:size(tetsToChange_1(:, 4), 1)
               if ismember(tetsToChange_1(numX, 4), X_IDs.bottomVerticesIds)
                   X(tetsToChange_1(numX, 4), :) = mean(X(X_IDs.bottomFaceIds(tetsToChange_1(numX, 1:3)), :));
               else
                   X(tetsToChange_1(numX, 4), :) = mean(X(X_IDs.topFaceIds(tetsToChange_1(numX, 1:3)), :));
               end
           end
%  
%            figure, tetramesh(tetsToChange_1, X);


            %Update Ys of intercalation cells
            %%% ONLY CHANGE 1 CELL (NOT DEBRIS)
%             changedYs = find(sum(ismember(tetrahedra, unique(tetsToChange_1)), 2)>2);
%             for numTetrahedron = changedYs'
%                 Y.DataRow(numTetrahedron, 1:2) = mean(X(tetrahedra(numTetrahedron, :), 1:2));
%             end
           
           correspondancePreviousNewTets = arrayfun(@(x) find(sum(ismember(tetrahedra_p(missingTets, :), tetrahedra(x, :)), 2)>2), tetrahedra_.n+1:tetrahedra_.n+size(newTetsModified, 1), 'UniformOutput', false);
           withoutCorrespondanceIds = missingTets(ismember(missingTets, missingTets(unique(vertcat(correspondancePreviousNewTets{:})))) == 0);
           numPrev = 1;
           for numTetrahedron = tetrahedra_.n+1:tetrahedra_.n+size(newTetsModified, 1)
               if isempty(correspondancePreviousNewTets{numPrev})
                   previousIdTet = withoutCorrespondanceIds(sum(ismember(tetrahedra_p(withoutCorrespondanceIds, :), tetrahedra(numTetrahedron, :)), 2) > 1);
               else
                   previousIdTet = missingTets(correspondancePreviousNewTets{numPrev});
               end
               oldTetCoords = mean(Yp.DataRow(previousIdTet, :), 1);
               newCoords = mean(X(tetrahedra(numTetrahedron, :), 1:3));
               newCoords = mean([newCoords; oldTetCoords], 1);
               newCoords(3) = oldTetCoords(3);
               Y = Y.Add(newCoords);
               Y0 = Y0.Add(newCoords);
               Yn = Yn.Add(newCoords);
               numPrev = numPrev + 1;
           end
           tetrahedra_ = tetrahedra_.Add(newTetsModified);
           
           %% 
           tetrahedra_ = tetrahedra_.RemoveCompletely(missingTets);
           Y = Y.RemoveCompletely(missingTets);
           Y0 = Y0.RemoveCompletely(missingTets);
           Yn = Yn.RemoveCompletely(missingTets);
          
          %Remove faces belonging to the cells in the intercalation
           faceToRemove = find(any(ismember(Cell.AllFaces.Nodes, unique(newTetsModified)), 2));
%            Cell.AllFaces=Cell.AllFaces.RemoveCompletely(faceToRemove);
%            SCn=SCn.RemoveCompletely(faceToRemove);
%            Cell.FaceCentres=Cell.FaceCentres.RemoveCompletely(faceToRemove);
           
           for numCellToChage = nodesToChange'
               currentFaces = Cell.Faces{numCellToChage};
               facesToRemove = ismember(currentFaces.FaceCentresID, faceToRemove);
               currentFaces.FaceCentresID(facesToRemove) = [];
               currentFaces.Tris(facesToRemove) = [];
               currentFaces.Vertices(facesToRemove) = [];
               currentFaces.nFaces = currentFaces.nFaces - sum(facesToRemove);
               Cell.Faces{numCellToChage} = currentFaces;
               Cell.cNodes{numCellToChage}(facesToRemove) = [];
           end

           %% Rebuild cells
           Cell.AssembleNodes=Cell.Int;
           [Cell, nC, SCn, flag]=ReBuildCells(Cell, tetrahedra_, Y, X, SCn);

           if ~flag
               Cell.AllFaces=Cell.AllFaces.CheckInteriorFaces(XgID);
               [Cell]=ComputeCellVolume(Cell,Y);
               Cell = Cell.computeEdgeLengths(Y);
               for jj=1:Cell.n
                   Cell.SAreaTrin{jj}=Cell.SAreaTri{jj};
                   Cell.EdgeLengthsn{jj}=Cell.EdgeLengths{jj};
               end
               Set.NumMainV=Y.n;
               Set.NumAuxV=Cell.FaceCentres.n;
               Set.NumCellCentroid = Cell.n;
               Set.NumTotalV=Set.NumMainV + Set.NumAuxV + Set.NumCellCentroid;
               [Cn]=BuildCn(tetrahedra_.Data);
               [Cell,Y]=CheckOrderingOfTriangulaiton(Cell,Y,Set);
           else
               return
           end
           
           if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.DataRow(1:tetrahedra_.n, :),Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end
           
           %% Solve modelling step with only those vertices
           [Dofs] = GetDOFs(Y, Cell, Set, isempty(Set.InputSegmentedImage) == 0);
           %changedYs = find(sum(ismember(tetrahedra_.DataRow, unique(tetsToChange_1)), 2)>2);
           [Dofs] = updateRemodelingDOFs(Dofs, Y.n - length(newTetsModified):Y.n, nC, Y);
           
           Cell.RemodelledVertices=[find(any(ismember(tetrahedra_.DataRow, nodesToChange), 2))', nC+Y.n];

           [Cell,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);
           
       end

       Cell.AssembleAll=true;
    end
end

