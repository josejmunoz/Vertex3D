function [Tnew, Ynew] = ConnectTetrahedra(Geo, nodeToRemove, nodesToChange, oldTets, mainNodes, Set, flipName, cellNodeLoosing)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

Tnew = [];
allTs = vertcat(Geo.Cells.T);
if isequal(Set.InputGeo, 'Voronoi')
    if length(cellNodeLoosing) == 1
        nodesLoosing = [cellNodeLoosing, mainNodes(vertcat(Geo.Cells(mainNodes).AliveStatus) == 0)];
    end
    
    mainNodesToConnect = setdiff(mainNodes, nodesLoosing);
    
    if length(mainNodes) >= 4
        nodesToConnect = unique([unique(allTs(sum(ismember(allTs, nodesLoosing), 2)> 1, :)); nodesToChange]);
        
        nodesConnectedToMainNodes = unique([getNodeNeighbours(Geo, mainNodesToConnect(1)); getNodeNeighbours(Geo, mainNodesToConnect(2))]);
        nodesConnectedToLoosingNodes = intersect(nodesToConnect, getNodeNeighbours(Geo, cellNodeLoosing));
        
        newCellBoundaryNode = setdiff(nodesConnectedToLoosingNodes, nodesConnectedToMainNodes);
        newCellBoundaryNode = newCellBoundaryNode(ismember(newCellBoundaryNode, Geo.XgID));
        
        if length(newCellBoundaryNode) > 1
            error('Need to check this!')
        end
        
        newCellBoundaryNode_Neighbours = getNodeNeighbours(Geo, newCellBoundaryNode);
        opposedNodesToConnect = setdiff(intersect(nodesConnectedToLoosingNodes, newCellBoundaryNode_Neighbours), [nodeToRemove, newCellBoundaryNode]);
        
        Tnew = [newCellBoundaryNode, nodeToRemove, mainNodesToConnect(1), intersect(getNodeNeighbours(Geo, mainNodesToConnect(1)), opposedNodesToConnect); ...
            newCellBoundaryNode, nodeToRemove, mainNodesToConnect(2), intersect(getNodeNeighbours(Geo, mainNodesToConnect(2)), opposedNodesToConnect)];
        
    else %% 3 mainNodes ('common')
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

