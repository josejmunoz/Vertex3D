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
        
        gNodeNeighbours = {};
        for numRow = 1:size(segmentFeatures, 1)
            gNodeNeighbours{numRow} = getNodeNeighbours(Geo, segmentFeatures{numRow, 2});
        end
        gNodes_NeighboursShared = unique(vertcat(gNodeNeighbours{:}));
        cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
        if sum([Geo.Cells(cellNodesShared).AliveStatus]) > 2 
            for numPair = 1:size(segmentFeatures, 1)
                if numPair == size(segmentFeatures, 1)
                    Set.NeedToConverge = 1;
                else
                    Set.NeedToConverge = 0;
                end
                
                hasConverged(numPair) = 0;

                cellNode = segmentFeatures{numPair, 1};
                ghostNode = segmentFeatures{numPair, 2};
                cellToIntercalateWith = segmentFeatures{numPair, 3};
                nodesShared = segmentFeatures{numPair, 5};
                nodesShared = nodesShared{1};
                faceGlobalId = segmentFeatures{numPair, 6};

                cellNodesShared = nodesShared(~ismember(nodesShared, Geo.XgID));
                ghostNodesShared = nodesShared(ismember(nodesShared, Geo.XgID));

                cellNodes = union(cellNodesShared, cellNode);

                %if ~all(ghostNodes) &&
                % If the shared nodes are all ghost nodes, we won't remodel

                %%if sum([Geo.Cells(cellNodes).AliveStatus]) >= 2 %&& ~any(ismember(faceGlobalId, newYgIds))
                nodesPair = [cellNode ghostNode];


                [valenceSegment, oldTets] = edgeValence(Geo, nodesPair);
                
                gNodes_NeighboursShared = unique(oldTets);
                cellNodesShared = gNodes_NeighboursShared(~ismember(gNodes_NeighboursShared, Geo.XgID));
                if sum([Geo.Cells(cellNodesShared).AliveStatus]) < 2
                    disp('hey')
                end

                %% Intercalation
                switch valenceSegment
                    case 2 %??
                        disp('error: valence tet 2')
                        sprintf('%s error: valence tet 2\n', Geo.log)
                        %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 3
                        disp('error: valence tet 3')
                        sprintf('%s error: valence tet 3\n', Geo.log)
                        %[Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip32(numFace, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 4
                        [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip4N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 5
                        [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip5N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    case 6
                        [Geo_0, Geo_n, Geo, Dofs, Set, newYgIds, hasConverged(numPair)] = Flip6N(nodesPair, oldTets, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds);
                    otherwise
                        disp('error: valence number greater than expected')
                        sprintf('%s error: valence number greater than expected\n', Geo.log)
                end
            end

            if any(~hasConverged)
                % Go back to initial state
                Geo_backup.log = Geo.log;
                Geo   = Geo_backup;
                Geo_n = Geo_n_backup;
                Dofs = Dofs_backup;
                Geo_0 = Geo_0_backup;
            end
        end
        
        checkedYgIds(end+1:end+size(segmentFeatures, 1), :) = [segmentFeatures{:, 1}, segmentFeatures{:, 2}];
        
        [segmentFeatures_all] = GetTrisToRemodelOrdered(Geo, Set);
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

