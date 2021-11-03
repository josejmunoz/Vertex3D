function [Cell, Y, tetrahedra] = simpleRemodelling(Cell, Y0, Yn, Y, CellInput, tetrahedra_, Tetrahedra_weights, X, X_IDs, SCn, XgID, Cn, verticesInfo, neighboursNetwork, Set)
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
           
           newConnectedNodes = nodesToChange(ismember(nodesToChange, currentIntercalation) == 0);
           
           % Reconstruct Tets
%            index = find(tetsToChangeAll_IDs)';
%            figure, tetramesh(tetrahedra(index, :), X);
           %% Option 1: New triangles
           % Change old connection for new connections
           trianglesConnectivity_all = verticesInfo.connectedCells;
           trianglesConnectivity_all(all(ismember(trianglesConnectivity_all, nodesToChange), 2), :) = sort([repmat(newConnectedNodes', 2, 1), sort(currentIntercalation)'], 2);
           neighboursNetwork(all(ismember(neighboursNetwork, currentIntercalation), 2), :) = [newConnectedNodes(1) newConnectedNodes(2); newConnectedNodes(2) newConnectedNodes(1)];
           
           % Remove vertices from edges
           vertices2DToChange = find(all(ismember(trianglesConnectivity_all, nodesToChange), 2));
           verticesConnectedToCells = currentIntercalation'; %sort(currentIntercalation)'; 
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
           changedTets  = ismember(newTets, tetrahedra, 'rows') == 0;
           tetrahedra = newTets;
           tetrahedra_ = DynamicArray(ceil(size(newTets,1)*1.5),size(newTets,2));
           tetrahedra_=tetrahedra_.Add(newTets);
           
%            index = find(tetsToChangeAll_IDs)';
%            figure, tetramesh(tetrahedra(index, :), X);
           
%            % New nodes
%            tetsToChange_1 = tetrahedra(sum(ismember(tetrahedra, nodesToChange), 2) >= 3 , :);
%            figure, tetramesh(tetsToChange_1, X);
%            for numX = 1:size(tetsToChange_1(:, 4), 1)
%                if ismember(tetsToChange_1(numX, 4), X_IDs.bottomVerticesIds)
%                    X(tetsToChange_1(numX, 4), :) = mean(X(X_IDs.bottomFaceIds(tetsToChange_1(numX, 1:3)), :));
%                else
%                    X(tetsToChange_1(numX, 4), :) = mean(X(X_IDs.topFaceIds(tetsToChange_1(numX, 1:3)), :));
%                end
%            end
 
%            figure, tetramesh(tetsToChange_1, X);
           
           % New Y_s %% TODO: HERE I NEED TO CHANGE THE EXACT TETS TO MAKE IT WORK
           Yp = Y;
           for numTetrahedron = 1:size(tetrahedra, 1) %find(changedTets)' % 
               newCoords = mean(X(tetrahedra(numTetrahedron, :), 1:2));
               Y.DataRow(numTetrahedron, 1:2) = newCoords;
           end
           
%           for numTetrahedron = 1:size(tetrahedra, 1)
%                if any(ismember(tetrahedra(numTetrahedron, :), verticesToChange))==0
%                    Y.DataRow(numTetrahedron, 1:2) = Yp.DataRow(numTetrahedron, :);
%                end
%            end
           
           %% Remove faces belonging to the cells in the intercalation
           faceToRemove = find(all(ismember(Cell.AllFaces.Nodes, nodesToChange), 2));
           Cell.AllFaces=Cell.AllFaces.Remove(faceToRemove);
           SCn=SCn.Remove(faceToRemove);
           Cell.FaceCentres=Cell.FaceCentres.Remove(faceToRemove);
           
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
           end
           
           %% Solve modelling step with only those vertices
           [Dofs] = GetDOFs(Y, Cell, Set, isempty(Set.InputSegmentedImage) == 0);
           [Dofs] = updateRemodelingDOFs(Dofs, 1:size(tetrahedra, 1), nC, Y);
           
           Cell.RemodelledVertices=[find(tetsToChangeAll_IDs)', nC+Y.n];
           if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.Data,Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end

           [Cell,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);
       end
    end
end

