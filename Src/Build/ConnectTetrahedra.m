function [Tnew, Ynew] = ConnectTetrahedra(Geo, nodeToRemove, nodesToChange, oldTets, mainNodes, Set, flipName, cellNodeLoosing)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

Tnew = [];
if isequal(Set.InputGeo, 'Voronoi')
    [closestDistance, closestID] = pdist2(vertcat(Geo.Cells(nodesToChange).X), Geo.Cells(nodeToRemove).X, 'euclidean', 'Smallest', 1);
    nodesToCombine = [nodesToChange(closestID), nodeToRemove];
    oldYs = cellfun(@(x) GetYFromTet(Geo, x), num2cell(oldTets, 2), 'UniformOutput', false);
    oldYs = vertcat(oldYs{:});
    if length(mainNodes) >= 4
        tetsToChange = oldTets(sum(ismember(oldTets, [nodesToChange(closestID) nodeToRemove]), 2) > 1, :);
        Ynew = oldYs(sum(ismember(oldTets, [nodesToChange(closestID) nodeToRemove]), 2) > 1, :);
        mainNodes(ismember(mainNodes, cellNodeLoosing)) = [];
        mainNodes(vertcat(Geo.Cells(mainNodes).AliveStatus) == 0) = [];
        if size(tetsToChange, 1) == 2
            tetXs = zeros(size(tetsToChange, 1), 3);
            for numTet = 1:size(tetsToChange, 1)
                tetXs(numTet, :) = mean(vertcat(Geo.Cells(tetsToChange(numTet, :)).X));
            end
            
            for nodeToConnect = mainNodes'
                nodeX = vertcat(Geo.Cells(nodeToConnect).X);
                [closestDistance, closestID] = pdist2(tetXs(:, 1:2), nodeX(1:2), 'euclidean', 'Smallest', 1);
                
                newTet = tetsToChange(closestID, :);
                newTet(newTet == cellNodeLoosing) = nodeToConnect;
                Tnew = vertcat(Tnew, newTet);
                tetsToChange(closestID, :) = [];
                tetXs(closestID, :) = [];
            end
        else
            error('WARNINGGGG check connect tetrahedra!!!!')
        end
    else %% 3 mainNodes ('common')
        [~, Tnew, Ynew, removedTets, replacedTets] = CombineTwoGhostNodes(Geo, Set, nodesToCombine, oldTets, oldYs);
    end
    
    Tnew
else
    
    if length(nodesToChange) > 4
        Tnew = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X), {'Qv', 'Q7'}));
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

