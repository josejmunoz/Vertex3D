function [Tnew] = ConnectTetrahedra(Geo, nodesToChange, oldTets, mainNodes, flipName)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

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

