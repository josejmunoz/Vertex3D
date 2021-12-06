function [Cell,Y0, Y,Yn,SCn,tetrahedra_,X,Dofs,Cn, Tetrahedra_weights, Set, verticesInfo, neighboursNetwork] = simpleRemodelling(Cell, Y0, Yn, Y, CellInput, tetrahedra_, X, X_IDs, SCn, XgID, Cn, verticesInfo, neighboursNetwork, Dofs, Tetrahedra_weights, Set)
%SIMPLEREMODELLING Summary of this function goes here
%   Detailed explanation goes here

    tetrahedra = tetrahedra_.DataRow;
    allVerticesRemodelled = [];
    allFacesRemodelled = [];

    %% Identify if any edge is shorter than it should
    [remodellingCells, apicalLengths,basalLengths] = identifyEdgesToIntercalate(Cell, tetrahedra, Set);
    
    while isempty(remodellingCells) == 0
       
       [smallestEdge, idEdges] = max(remodellingCells(:, 3:4), [], 1);
       [~, idSmallestEdge] = max(smallestEdge);
       remodellingCells = remodellingCells(idEdges(idSmallestEdge), :);
       
       %% Intercalate whole cell in 3D: both cells do not share a face anymore 
       for numIntercalation = 1:size(remodellingCells, 1)
           tetrahedra_p = tetrahedra;
           Cellp = Cell; Setp = Set; CellInputp = CellInput;
           Yp = Y; SCnp = SCn; Y0p = Y0; Ynp = Yn;
           currentIntercalation = remodellingCells(numIntercalation, 1:2); % Face to remove
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
           verticesInfo.connectedCells = trianglesConnectivity_all;
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
           
           [Twg_bottom] = createTetrahedra(trianglesConnectivity_all, neighboursNetwork, verticesInfo.edges, Cell.Int', X_IDs.bottomFaceIds, X_IDs.bottomVerticesIds);
           [Twg_top] = createTetrahedra(trianglesConnectivity_all, neighboursNetwork, verticesInfo.edges, Cell.Int', X_IDs.topFaceIds, X_IDs.topVerticesIds);
           newTets = vertcat(Twg_top, Twg_bottom);
           newTets(all(ismember(newTets,XgID),2),:) = [];
           [changedTets, Locb] = ismember(tetrahedra, newTets, 'rows');
           [changedTets_2, Locb] = ismember(newTets, tetrahedra, 'rows');
           tetIds = 1:tetrahedra_.n;
           missingTets = tetIds(ismember(tetIds, Locb) == 0);
           newTetsModified = newTets(changedTets_2==0, :);
           
           tetrahedra(missingTets, :) = repmat([0 0 0 0], length(missingTets), 1);
           tetrahedra_.DataRow(missingTets, :) = repmat([0 0 0 0], length(missingTets), 1);
           Y.DataRow(missingTets, :) = repmat([-100 -100 -100], length(missingTets), 1); %% these vertices don't affect the mechanics
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
           
           correspondancePreviousNewTets = arrayfun(@(x) find(sum(ismember(tetrahedra_p(missingTets, :), tetrahedra(x, :)), 2)>2), tetrahedra_.n+1:tetrahedra_.n+size(newTetsModified, 1), 'UniformOutput', false);
           withoutCorrespondanceIds = missingTets(ismember(missingTets, missingTets(unique(vertcat(correspondancePreviousNewTets{:})))) == 0);
           numPrev = 1;
           newVerticesIDs = tetrahedra_.n+1:tetrahedra_.n+size(newTetsModified, 1);
           for numTetrahedron = newVerticesIDs
               if isempty(correspondancePreviousNewTets{numPrev})
                   previousIdTet = withoutCorrespondanceIds(sum(ismember(tetrahedra_p(withoutCorrespondanceIds, :), tetrahedra(numTetrahedron, :)), 2) > 1);
               else
                   previousIdTet = missingTets(correspondancePreviousNewTets{numPrev});
                   
               end
               oldTetCoords = mean(Yp.DataRow(previousIdTet, :), 1);
               newCoords = mean(X(tetrahedra(numTetrahedron, :), 1:3));
               newCoords = newCoords.*0.8 + oldTetCoords .*0.2;
               newCoords(3) = oldTetCoords(3);
               Y = Y.Add(newCoords);
               Y0 = Y0.Add(newCoords);
               Yn = Yn.Add(newCoords);
               Tetrahedra_weights(end+1, :) = mean(Tetrahedra_weights(previousIdTet, :), 1);
               numPrev = numPrev + 1;
           end
           tetrahedra_ = tetrahedra_.Add(newTetsModified);
           Tetrahedra_weights(missingTets, :) = repmat([0 0 0], length(missingTets), 1);

           %% Rebuild cells
           Cell.AssembleNodes=nodesToChange; %Cell.Int;
           [Cell, newFaces, SCn, flag]=ReBuildCells(Cell, tetrahedra_, Y, X, SCn);

           if ~flag
               allFaces = [Cell.Faces{:}];
               usedIDFaces = unique(vertcat(allFaces.FaceCentresID));
               unusedIDFaces = setdiff(1:max(usedIDFaces), usedIDFaces);
               Cell.FaceCentres.DataRow(unusedIDFaces, :) = repmat([-100 -100 -100], length(unusedIDFaces), 1);
               Cell.FaceCentres0.DataRow(unusedIDFaces, :) =  Cell.FaceCentres.DataRow(unusedIDFaces, :);
               SCn.DataRow(unusedIDFaces, :) = Cell.FaceCentres.DataRow(unusedIDFaces, :);
               
               Cell.AllFaces=Cell.AllFaces.CheckInteriorFaces(XgID);
               [Cell]=ComputeCellVolume(Cell,Y);
               Cell = Cell.computeEdgeLengths(Y);
               [Cell] = Cell.computeEdgeLocation(Y);
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
               %% ERROR
               error('Error rebuilding cells at remodelling');
           end

           %if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.DataRow(1:tetrahedra_.n, :),Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end
           
           %% Solve modelling step with only those vertices
           allVerticesRemodelled{end+1} = nodesToChange;
           allFacesRemodelled{end+1} = newFaces;
           
           [remodellingCells, apicalLengths,basalLengths] = identifyEdgesToIntercalate(Cell, tetrahedra, Set);
           
           if isempty(remodellingCells) == 0
               invalidRemodellingCells = cellfun(@(x) find(sum(ismember(remodellingCells(:, 1:2), x), 2) == 2)', allVerticesRemodelled, 'UniformOutput', false);
               invalidRemodellingCells = unique([invalidRemodellingCells{:}]);
               remodellingCells(invalidRemodellingCells, :) = [];
           end
           if isempty(remodellingCells)
               if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.DataRow(1:tetrahedra_.n, :),Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr+1,Set); end
               nodesToChange = unique([allVerticesRemodelled{:}]);
               newFaces = unique([allFacesRemodelled{:}]);
               Cell.RemodelledVertices=find(sum(ismember(tetrahedra_.Data, nodesToChange), 2) > 0);
               changedFaces = find(any(ismember(Cell.AllFaces.Nodes, nodesToChange), 2));
               remodelledFaces = setdiff([newFaces changedFaces'], unusedIDFaces);

               [Dofs] = GetDOFs(Y, Cell, Set, isempty(Set.InputSegmentedImage) == 0, tetrahedra_.DataRow);
               [Dofs] = updateRemodelingDOFs(Dofs, Cell.RemodelledVertices, remodelledFaces, Y);

               Y0.DataRow(Cell.RemodelledVertices, :) = Y.DataRow(Cell.RemodelledVertices, :);
               Cell.FaceCentres0.DataRow(remodelledFaces, :) = Cell.FaceCentres.DataRow(remodelledFaces, :);
               
               maxSteps = 10;
               Set.nu = Setp.nu * 100;
               Set.nu0 = Setp.nu * 20;
               for numStep = 1:maxSteps
                   numStep

                   SCn = Cell.FaceCentres;
                   Cell.Centre_n = Cell.Centre;                       
                   Yn=Y;
                   
                   [Set] = updateMechanicalParams(Set, Setp, maxSteps, numStep);
                   Cell.Vol0 = Cell.Vol .* ((maxSteps - numStep + 1)/maxSteps) + Cellp.Vol0 .* ((numStep - 1)/maxSteps);
                   Cell.Centre0 = Cell.Centre .* ((maxSteps - numStep + 1)/maxSteps) + Cellp.Centre0 .* ((numStep - 1)/maxSteps);
                   [Cell,Y,Yn,SCn,X,Dofs,Set,~,DidNotConverge]=SolveRemodelingStep(Cell,Y0,Y,X,Dofs,Set,Yn,SCn,CellInput);

                   if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.DataRow(1:tetrahedra_.n, :),Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr+numStep+1,Set); end

                   Set.MaxIter=Set.MaxIter0;
                   Set.ApplyBC=true;
                   if numStep >= 3
                       [Dofs] = GetDOFs(Y, Cell, Set, isempty(Set.InputSegmentedImage) == 0, tetrahedra_.DataRow);
                       Dofs.Remodel = Dofs.FreeDofs;
                       Set.nu = Setp.nu * 1000;
                       Set.nu0 = Setp.nu * ((20*maxSteps) - (20*(numStep)));
                   end
               end
               Set = Setp;
               Yn=Y;
               SCn=Cell.FaceCentres;
               Set.NumMainV=Y.n;
               Set.NumAuxV=Cell.FaceCentres.n;
               Set.NumCellCentroid = Cell.n;
               Set.NumTotalV=Set.NumMainV + Set.NumAuxV + Set.NumCellCentroid;
               Y0.DataRow(Cell.RemodelledVertices, :) = Y.DataRow(Cell.RemodelledVertices, :);
               Cell.FaceCentres0.DataRow(remodelledFaces, :) = Cell.FaceCentres.DataRow(remodelledFaces, :);
           end
       end
       %Set.nu = Set.nu * 100;
       Set.ReModel=false;
       Cell.AssembleAll=true;        
       %if Set.VTK, PostProcessingVTK(X,Y,tetrahedra_.DataRow(1:tetrahedra_.n, :),Cn,Cell,strcat(Set.OutputFolder,Esc,'ResultVTK'),Set.iIncr,Set); end    
    end
end

