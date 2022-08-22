function [Tnew, Ynew] = ConnectTetrahedra(Geo, Geo_n, nodesToChange, oldTets, mainNode, Set, flipName)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

nodeToRemove = setdiff(unique(oldTets), nodesToChange);

Tnew = [];
if isequal(Set.InputGeo, 'Voronoi')
    [closestDistance, closestID] = pdist2(vertcat(Geo.Cells(nodesToChange).X), Geo.Cells(nodeToRemove).X, 'euclidean', 'Smallest', 1);
    nodesToCombine = [nodesToChange(closestID), nodeToRemove];
    oldYs = cellfun(@(x) GetYFromTet(Geo, x), num2cell(oldTets, 2), 'UniformOutput', false);
    oldYs = vertcat(oldYs{:});
    [~, Tnew, Ynew, removedTets, replacedTets] = CombineTwoGhostNodes(Geo, Set, nodesToCombine, oldTets, oldYs);
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
        Tnew = horzcat(ones(size(DT.ConnectivityList, 1), 1) * mainNode, nodesToChange(DT.ConnectivityList));
    end
end
end

