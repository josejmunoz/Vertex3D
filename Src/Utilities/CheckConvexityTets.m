function [convexTet] = CheckConvexityTets(nodeTets, newTets, Geo, nodeItShouldBeConnected)
%CHECKCONVEXITYTETS Summary of this function goes here
%   Detailed explanation goes here

%Combine the two new tets created
uniqueNodes = unique([nodeTets; newTets]);
%Obtain the CONCAVE shape connecting the tets (this is to allow concave
%shapes between different tets)
shapeBoundary = boundary(vertcat(Geo.Cells(uniqueNodes).X), 1);
uniqueConnectedNodes = uniqueNodes(shapeBoundary);
%Obtain how they are connected
neighbours = unique(uniqueConnectedNodes(any(ismember(uniqueConnectedNodes, nodeItShouldBeConnected), 2), :));


%Visualize if it is correct
% allNodes = vertcat(Geo.Cells.X);
% figure, trisurf(uniqueConnectedNodes, allNodes(:, 1), allNodes(:, 2), allNodes(:, 3))
% text(allNodes(uniqueNodes, 1), allNodes(uniqueNodes, 2), allNodes(uniqueNodes, 3), cellfun(@num2str, num2cell(uniqueNodes), 'UniformOutput', false),'VerticalAlignment','bottom','HorizontalAlignment','right')

%It would be convex if the connections correspond to their assigned node
%(nodeTets)
convexTet = all(ismember(nodeTets, neighbours));
end

