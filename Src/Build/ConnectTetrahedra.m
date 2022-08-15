function [Tnew] = ConnectTetrahedra(Geo, nodesToChange, oldTets, mainNode, flipName)
%CONNECTTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

if length(nodesToChange) > 4
    Tnew = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X)));
else
    Tnew = nodesToChange';
end

% Remove tets with all Ghost Nodes
Tnew(all(ismember(Tnew, Geo.XgID), 2), :) = [];

%% Check if everything is correct and try to correct otherwise
[overlappingTets, correctedTets] = CheckOverlappingTets(oldTets, Tnew, Geo, flipName);

% if ~isempty(correctedTets)
%     Tnew = correctedTets;
% end

if length(nodesToChange) > 4 && overlappingTets
    %% NEED TO DO THIS INSTEAD: https://meshlogic.github.io/posts/jupyter/curve-fitting/fitting-a-circle-to-cluster-of-3d-points/
    nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus})) = [];
    [~,score] = pca(vertcat(Geo.Cells(nodesToChange).X));
    DT = delaunayTriangulation(score(:, 1:2));
    Tnew = horzcat(ones(size(DT.ConnectivityList, 1), 1) * mainNode, nodesToChange(DT.ConnectivityList));
end
end

