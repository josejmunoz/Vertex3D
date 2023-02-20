function [Geo_0, Geo_n, Geo, Dofs, Set] = Remodeling(Geo_0, Geo_n, Geo, Dofs, Set)

    
    Geo.AssemblegIds = [];
    newYgIds = [];
    checkedYgIds = [];
    
    [segmentFeatures_all] = GetTrisToRemodelOrdered(Geo, Set);
    
    %% loop ENERGY-dependant
    while ~isempty(segmentFeatures_all)
        Geo_backup = Geo; Geo_n_backup = Geo_n; Geo_0_backup = Geo_0; Dofs_backup = Dofs;
        
        segmentFeatures = segmentFeatures_all{1};
        [~, ids] = unique(segmentFeatures(:, 1:2), 'rows');
        segmentFeatures = segmentFeatures(ids, :);
        segmentFeatures = sortrows(segmentFeatures, 6);
        gNodeNeighbours = {};
        for numRow = 1:size(segmentFeatures, 1)
            gNodeNeighbours{numRow} = getNodeNeighbours(Geo, segmentFeatures{numRow, 2});
        end
        gNodes_NeighboursShared = unique(vertcat(gNodeNeighbours{:}));
        cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
        if sum([Geo.Cells(cellNodesShared).AliveStatus]) > 2 
            Set.NeedToConverge = 0;
            allTnew = [];
            initialNodeValence = arrayfun(@(x) sum((ismember(getNodeNeighbours(Geo, x), Geo.XgID))), [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]);
            for numPair = 1:size(segmentFeatures, 1)
                
                cellNode = segmentFeatures{numPair, 1};
                ghostNode = segmentFeatures{numPair, 2};
                cellToIntercalateWith = segmentFeatures{numPair, 3};
                cellToSplitFrom = segmentFeatures{numPair, 4};
                
                nodesToCombineLater = [];
                hasConverged(numPair) = 1;
                while hasConverged(numPair) == 1
                    hasConverged(numPair) = 0;

                    nodesPairs = [cellNode ghostNode];
                    
                    nodesToCombineLater(end+1) = ghostNode;
                    
                    for nodesPair = nodesPairs'
                    
                        [valenceSegment, oldTets, oldYs] = edgeValence(Geo, nodesPair);

                        %% Intercalation
                        switch valenceSegment
                            case 0
                                break;
                            case 2 %??
                                disp('error: valence tet 2')
                                sprintf('%s error: valence tet 2\n', Geo.log)
                                %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                                break;
                            case 3
                                disp('error: valence tet 3')
                                sprintf('%s error: valence tet 3\n', Geo.log)
                                break;
                                %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip32(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                            case 4
                                [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip4N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                            case 5
                                [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip5N(nodesPair, oldTets, oldYs, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                            case 6
                                [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip6N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                            case 7
                                [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair), Tnew] = Flip7N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                            otherwise
                                disp('error: valence number greater than expected')
                                sprintf('%s error: valence number greater than expected\n', Geo.log);
                                break;
                        end
                    
                        allTnew = vertcat(allTnew, Tnew);
                    end

                    sharedNodesStill = getNodeNeighboursPerDomain(Geo, cellNode, ghostNode, cellToSplitFrom);

                    if any(ismember(sharedNodesStill, Geo.XgID))
                        sharedNodesStill_g = sharedNodesStill(ismember(sharedNodesStill, Geo.XgID));
                        ghostNode = sharedNodesStill_g(1);
                    else
                        PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
                        nonDeadCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
                        
                        remodellingNodeValence = arrayfun(@(x) sum((ismember(getNodeNeighbours(Geo, x), Geo.XgID))), [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]);
                        nodesToAddOrRemove = remodellingNodeValence - initialNodeValence;
                        
                        for nodeToChangeValence = find(nodesToAddOrRemove ~= 0)
                            for numTime = 1:abs(nodesToAddOrRemove(nodeToChangeValence))
                                if nodesToAddOrRemove(nodeToChangeValence) > 0
                                    [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = FlipN0(Geo, Geo_n, Geo_0, Dofs, newYgIds, nodeToRemove, nodeToKeep, Set);
                                else
                                    closerNode = getNodeNeighboursPerDomain(Geo, nodesToCombineLater(numTime), nodesToCombineLater(numTime));
                                    closerNode_g = closerNode(ismember(closerNode, Geo.XgID));
                                    nodesOfCellToAdd = getNodeNeighboursPerDomain(Geo, nonDeadCells(nodeToChangeValence), nodesToCombineLater(numTime));
                                    nodesOfCellToAdd_g = nodesOfCellToAdd(ismember(nodesOfCellToAdd, Geo.XgID));
                                    
                                    nodeCloserToCell = intersect(closerNode_g, nodesOfCellToAdd_g);
                                    if isempty(nodeCloserToCell)
                                        error('caca');
                                    end
                                    
                                    nodesToPick = getNodeNeighboursPerDomain(Geo, nodeCloserToCell(1), nodesToCombineLater(numTime), nonDeadCells(nodeToChangeValence));
                                    nodesToPick_g = nodesToPick(ismember(nodesToPick, Geo.XgID));
                                    
                                    edgeToBeAddedNode = [nodeCloserToCell(1), nodesToPick_g(1)];
                                    newNodes = mean(vertcat(Geo.Cells(edgeToBeAddedNode).X));
                                    [~, oldTets, ~] = edgeValence(Geo, edgeToBeAddedNode);
                                    oldTets = oldTets(any(ismember(oldTets, nonDeadCells(nodeToChangeValence)), 2), :);
                                    
                                    surroundingNodes = unique(oldTets(:));
                                    surroundingNodes = setdiff(surroundingNodes, setdiff(nonDeadCells, nonDeadCells(nodeToChangeValence)));
                                    [Geo, Geo_n, Geo_0, Dofs, Set, newYgIds, hasConverged] = FlipAddNodes(surroundingNodes, oldTets, newNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                                end
                            end
                        end
                        PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2);
                        
                        break;
                    end
                end
            end 

            %% Vertices connecting the two intercalating cells should be closer
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
            
            closeToNewPoint = 0.4;
            
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
            
            Geo   = Rebuild(Geo, Set);
            Geo   = BuildGlobalIds(Geo);
            Geo   = UpdateMeasures(Geo);
            Geo_n = Geo;
            Geo_0 = Rebuild(Geo_0, Set);
            Geo_0 = BuildGlobalIds(Geo_0);    
            PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
            
            %% 
            Dofs = GetDOFs(Geo, Set);
            [Dofs, Geo]  = GetRemodelDOFs(allTnew, Dofs, Geo);
            
            [g, K, E, ~, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
            dyr=norm(-K(Dofs.Remodel,Dofs.Remodel)\g(Dofs.Remodel))
            gr=norm(g(Dofs.Remodel))
            [g, K, E, ~, Energies_backup] = KgGlobal(Geo_0_backup, Geo_n_backup, Geo_backup, Set);
            dyr=norm(-K(Dofs.Remodel,Dofs.Remodel)\g(Dofs.Remodel))
            gr=norm(g(Dofs.Remodel))
            Energies_backup
            Energies
            [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
            if DidNotConverge
                % Go back to initial state
                Geo_backup.log = Geo.log;
                Geo   = Geo_backup;
                Geo_n = Geo_n_backup;
                Dofs = Dofs_backup;
                Geo_0 = Geo_0_backup;
                Geo.log = sprintf('%s =>> %s-Flip rejected: did not converge\n', Geo.log, flipName);
                return
            end

            newYgIds = unique([newYgIds; Geo.AssemblegIds]);
            Geo   = UpdateMeasures(Geo);

            PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)

            hasConverged = 1;
        end
        
        checkedYgIds(end+1:end+size(segmentFeatures, 1), :) = [segmentFeatures{:, 1}, segmentFeatures{:, 2}];
        
        %[segmentFeatures_all] = GetTrisToRemodelOrdered(Geo, Set);
        rowsToRemove = [];
        if ~isempty(segmentFeatures_all)
            for numRow = 1:length(segmentFeatures_all)
                cSegFea = segmentFeatures_all{numRow};
                if all(ismember([cSegFea{:, 1:2}], checkedYgIds, 'rows'))
                    rowsToRemove(end+1) = numRow;
                end
            end
        end
        segmentFeatures_all(rowsToRemove) = [];
    end
    
    [g, K, E, Geo, Energies] = KgGlobal(Geo_0, Geo_n, Geo, Set);
end

