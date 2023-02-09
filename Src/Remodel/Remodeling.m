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
        segmentFeatures = sortrows(sortrows(segmentFeatures, 3, 'descend'), 4);
        gNodeNeighbours = {};
        for numRow = 1:size(segmentFeatures, 1)
            gNodeNeighbours{numRow} = getNodeNeighbours(Geo, segmentFeatures{numRow, 2});
        end
        gNodes_NeighboursShared = unique(vertcat(gNodeNeighbours{:}));
        cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
        if sum([Geo.Cells(cellNodesShared).AliveStatus]) > 2 
            Set.NeedToConverge = 0;
            allTnew = [];
            for numPair = 1:size(segmentFeatures, 1)
                
                cellNode = segmentFeatures{numPair, 1};
                ghostNode = segmentFeatures{numPair, 2};
                cellToIntercalateWith = segmentFeatures{numPair, 3};
                cellToSplitFrom = segmentFeatures{numPair, 4};
                
                hasConverged(numPair) = 1;
                while hasConverged(numPair) == 1
                    hasConverged(numPair) = 0;

                    %if ~all(ghostNodes) &&
                    % If the shared nodes are all ghost nodes, we won't remodel

                    %%if sum([Geo.Cells(cellNodes).AliveStatus]) >= 2 %&& ~any(ismember(faceGlobalId, newYgIds))
                    nodesPair = [cellNode ghostNode];

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

                    sharedNodesStill = getNodeNeighboursPerDomain(Geo, cellNode, ghostNode, cellToSplitFrom);

                    if any(ismember(sharedNodesStill, Geo.XgID))
                        sharedNodesStill_g = sharedNodesStill(ismember(sharedNodesStill, Geo.XgID));
                        ghostNode = sharedNodesStill_g(1);
                    else
                        break;
                    end
                end
            end

            Geo   = Rebuild(Geo, Set);
            Geo   = BuildGlobalIds(Geo);
            Geo   = UpdateMeasures(Geo);
            Geo_n = Rebuild(Geo_n, Set);
            Geo_n = BuildGlobalIds(Geo_n);
            Geo_0 = Rebuild(Geo_0, Set);
            Geo_0 = BuildGlobalIds(Geo_0);   
            
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
            
            for tetToCheck = verticesToChange'
                for nodeInTet = tetToCheck'
                    if ~ismember(nodeInTet, Geo.XgID)
                        newPoint = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);

                        vectorRefNew = refPoint - newPoint;

                        Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = refPoint - vectorRefNew/2;
                        %Geo_n.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :) = Geo.Cells(nodeInTet).Y(ismember(sort(Geo.Cells(nodeInTet).T, 2), tetToCheck', 'rows'), :);
                    end
                end
            end
            
%             error('continue improving this!');
%             ghostNodesToUpdate = unique(verticesToChange(:));
%             ghostNodesToUpdate = ghostNodesToUpdate(ismember(ghostNodesToUpdate, Geo.XgID));
%             
%             for ghostNodeToUpdate = ghostNodesToUpdate'
%                 neighbours = unique(Geo.Cells(ghostNodeToUpdate).T);
%                 neighbours_Cells = neighbours(~ismember(neighbours, Geo.XgID));
% 
%                  newXValue = mean(vertcat(Geo.Cells(neighbours_Cells).X));
%                  Geo.Cells(ghostNodeToUpdate).X(1:2) = newXValue(1:2);
%             end
%             
%             Geo = BuildXFromY(Geo_n, Geo);
%             
%             %%% SHOULD I MOVE THE GHOST NODES ALSO??? FOR THE FACES TO BE
%             %%% WELL PLACED??
            Geo   = Rebuild(Geo, Set);
            Geo   = BuildGlobalIds(Geo);
            Geo   = UpdateMeasures(Geo);
            Geo_n = Rebuild(Geo_n, Set);
            Geo_n = BuildGlobalIds(Geo_n);
            Geo_n = UpdateMeasures(Geo_n);
            Geo_0 = Rebuild(Geo_0, Set);
            Geo_0 = BuildGlobalIds(Geo_0);    
%             %% Update Geo_0 to be reset the vertices that we have changed averaging with previous Geo_0 and current Geo
%             percentageGeo = 1 - Set.Reset_PercentageGeo0;
%             for c=1:Geo.nCells
%                 if ismember(c, Tnew) && ~isempty(Geo.Cells(c).AliveStatus) && Geo.Cells(c).AliveStatus == 1
%                     Geo_0.Cells(c).X = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).X + percentageGeo * Geo.Cells(c).X;
%                     Geo_0.Cells(c).Y = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).Y + percentageGeo * Geo.Cells(c).Y;
%                     
%                     for f=1:length(Geo.Cells(c).Faces)
%                         Geo_0.Cells(c).Faces(f).Centre = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).Faces(f).Centre + percentageGeo * Geo.Cells(c).Faces(f).Centre;
%                     end
%                 end
%             end
            
%             Geo_0 = UpdateMeasures(Geo_0);
            PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1);
            
            % Also the vertex middle Scutoid vertex
            
            
            %% 
            Dofs = GetDOFs(Geo, Set);
            [Dofs, Geo]  = GetRemodelDOFs(allTnew, Dofs, Geo);
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
            Geo_n = UpdateMeasures(Geo_n);

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

