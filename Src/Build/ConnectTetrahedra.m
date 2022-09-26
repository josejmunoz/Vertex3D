function [Tnew, Ynew, oldTets] = ConnectTetrahedra(Geo, nodeToRemove, nodesToChange, oldTets, mainNodes, Set, flipName, cellNodeLoosing)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

Tnew = [];
Ynew = [];

allTs = vertcat(Geo.Cells.T);
if isequal(Set.InputGeo, 'Voronoi')
    
    if length(mainNodes) == 4
        ghostNodeLoosing = mainNodes(vertcat(Geo.Cells(mainNodes).AliveStatus) == 0);

        if length(mainNodes(vertcat(Geo.Cells(mainNodes).AliveStatus) == 0)) == 1        
            if length(cellNodeLoosing) == 1
                nodesLoosing = [cellNodeLoosing, mainNodes(vertcat(Geo.Cells(mainNodes).AliveStatus) == 0)];
            end

            mainNodesToConnect = setdiff(mainNodes, nodesLoosing);
            nodesToConnect = unique([unique(allTs(sum(ismember(allTs, nodesLoosing), 2)> 1, :)); nodesToChange]);

            nodesConnectedToMainNodes = unique([getNodeNeighbours(Geo, mainNodesToConnect(1)); getNodeNeighbours(Geo, mainNodesToConnect(2))]);
            nodesConnectedToLoosingNodes = intersect(nodesToConnect, getNodeNeighbours(Geo, cellNodeLoosing));

            newCellBoundaryNode = setdiff(nodesConnectedToLoosingNodes, nodesConnectedToMainNodes);

            if ismember(nodeToRemove, Geo.XgTop)
                newCellBoundaryNode = newCellBoundaryNode(ismember(newCellBoundaryNode, Geo.XgTop));
            elseif ismember(nodeToRemove, Geo.XgBottom)
                newCellBoundaryNode = newCellBoundaryNode(ismember(newCellBoundaryNode, Geo.XgBottom));
            else
                newCellBoundaryNode = newCellBoundaryNode(ismember(newCellBoundaryNode, Geo.XgID));
            end
            
            if length(newCellBoundaryNode) > 1
                nodeNeighbours_Boundary = arrayfun(@(x) sum(~ismember(getNodeNeighbours(Geo, x), Geo.XgID)) == 1, newCellBoundaryNode);
                
                newCellBoundaryNode = newCellBoundaryNode(nodeNeighbours_Boundary);
                if length(newCellBoundaryNode) > 1
                    error('Need to check this!')
                end
            end
            
            connectedNodes = [nodeToRemove, newCellBoundaryNode];
            newCellBoundaryNode_Neighbours = getNodeNeighbours(Geo, newCellBoundaryNode);
            opposedNodesToConnect = setdiff(intersect(nodesConnectedToLoosingNodes, newCellBoundaryNode_Neighbours), [nodeToRemove, newCellBoundaryNode]);

            %% Connections #1: 1 mainNodes and 3 ghost node
            Tnew = [newCellBoundaryNode, nodeToRemove, mainNodesToConnect(1), intersect(getNodeNeighbours(Geo, mainNodesToConnect(1)), opposedNodesToConnect); ...
                newCellBoundaryNode, nodeToRemove, mainNodesToConnect(2), intersect(getNodeNeighbours(Geo, mainNodesToConnect(2)), opposedNodesToConnect)];

            %% Connections #2: 2 mainNodes and 2 ghost node
            Tnew(end+1, :) = [mainNodesToConnect', newCellBoundaryNode, nodeToRemove];
            Tnew(end+1, :) = [cellNodeLoosing, newCellBoundaryNode, mainNodesToConnect(1), intersect(getNodeNeighbours(Geo, mainNodesToConnect(1)), opposedNodesToConnect)];
            Tnew(end+1, :) = [cellNodeLoosing, newCellBoundaryNode, mainNodesToConnect(2), intersect(getNodeNeighbours(Geo, mainNodesToConnect(2)), opposedNodesToConnect)];

            %% Connections #3: 3 mainNodes and 1 ghost node
            connectedNode1_ghostNodeLoosing = intersect(getNodeNeighbours(Geo, ghostNodeLoosing), connectedNodes);
            connectedNode2_mainNodeLoosing = setdiff(connectedNodes, connectedNode1_ghostNodeLoosing);
            Tnew(end+1, :) = [mainNodesToConnect', ghostNodeLoosing, connectedNode1_ghostNodeLoosing];
            Tnew(end+1, :) = [mainNodesToConnect', cellNodeLoosing, connectedNode2_mainNodeLoosing];
            
            %% Connection #4: 4 mainNodes
            Tnew(end+1, :) = mainNodes;
            
            %% Update oldTets
            nodesChanged = unique(Tnew(:));
            oldTets = oldTets(sum(ismember(oldTets, nodesChanged), 2) > 3, :);
            
            %% Recalculate Ys
            allTs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).T);
            allYs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).Y);
            for numTet = 1:size(Tnew, 1)
                tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 2;
                
                mainNode_current = mainNodesToConnect(ismember(mainNodesToConnect, Tnew(numTet, :)));
                if any(tetsToUse)
                    contributionOldYs = 1;
                    Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * ComputeY(vertcat(Geo.Cells(Tnew(numTet, :)).X), Geo.Cells(mainNode_current(1)).X, length([Geo.Cells(Tnew(numTet, :)).AliveStatus]), Set);
                else
                    contributionOldYs = 0.9;
                    tetsToUse = sum(ismember(allTs, Tnew(numTet, :)), 2) > 1;
                    
                    if any(ismember(Tnew(numTet, :), Geo.XgTop))
                        tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
                    elseif any(ismember(Tnew(numTet, :), Geo.XgBottom))
                        tetsToUse = tetsToUse & any(ismember(allTs, Geo.XgTop), 2);
                    end
                    
                    if any(tetsToUse)
                        Ynew(end+1, :) = contributionOldYs * mean(vertcat(allYs(tetsToUse, :)), 1) + (1-contributionOldYs) * ComputeY(vertcat(Geo.Cells(Tnew(numTet, :)).X), Geo.Cells(mainNode_current(1)).X, length([Geo.Cells(Tnew(numTet, :)).AliveStatus]), Set);
                    else
                        error('Need to check this!');
                    end
                end
            end
        else
            %%error('Need to check this!');
        end
    else %% 3 mainNodes ('common')
        [~, closestID] = pdist2(vertcat(Geo.Cells(nodesToChange).X), Geo.Cells(nodeToRemove).X, 'euclidean', 'Smallest', 1);
        nodesToCombine = [nodesToChange(closestID), nodeToRemove];
        oldYs = cellfun(@(x) GetYFromTet(Geo, x), num2cell(oldTets, 2), 'UniformOutput', false);
        oldYs = vertcat(oldYs{:});
        [~, Tnew, Ynew, removedTets, replacedTets] = CombineTwoGhostNodes(Geo, Set, nodesToCombine, oldTets, oldYs);
    end
else
    
    if length(nodesToChange) > 4
        Tnew = nodesToConnect(delaunayn(vertcat(Geo.Cells(nodesToConnect).X), {'Qv', 'Q7'}));
    else
        Tnew = nodesToChange';
    end

    % Remove tets with all Ghost Nodes
    Tnew(all(ismember(Tnew, Geo.XgID), 2), :) = [];

    %% Check if everything is correct and try to correct otherwise
    [overlappingTets, correctedTets] = CheckOverlappingTets(oldTets, Tnew, Geo, flipName);

    if ~isempty(correctedTets)
        Tnew = correctedTets;
        [overlappingTets] = CheckOverlappingTets(oldTets, Tnew, Geo, flipName);
    end

    if length(nodesToChange) > 4 && overlappingTets && sum(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus})) == 1
        %% NEED TO DO THIS INSTEAD: https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
        nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus})) = [];
        [~,score] = pca(vertcat(Geo.Cells(nodesToChange).X));
        DT = delaunayTriangulation(score(:, 1:2));
        Tnew = horzcat(ones(size(DT.ConnectivityList, 1), 1) * mainNodes, nodesToChange(DT.ConnectivityList));
    end
end
end

