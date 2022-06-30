function [isConvex] = CheckConvexityTets(nodeTets, Tnew, Geo, nodeItShouldBeConnected)
%CHECKCONVEXITYTETS Summary of this function goes here
%   Detailed explanation goes here

isConvex = true;
%Get the centre of the tetrahedron Tnew and look for if it is inside any
%other tet.
for numTnew = 1:size(Tnew, 1)
    newTet = Tnew(numTnew, :);
    newTetCentroid = mean(vertcat(Geo.Cells(newTet).X), 1);
    
    tetXs = vertcat(Geo.Cells(nodeTets).X);
    tetShape = alphaShape(tetXs);
    
    % Checking if any point of the Xs are inside the tetrahedra
    if any(tetShape.inShape(newTetCentroid(1), newTetCentroid(2), newTetCentroid(3)))
%         figure, plot(tetShape)
%         text(tetXs(:, 1), tetXs(:, 2), tetXs(:, 3), cellfun(@num2str, num2cell(nodeTets), 'UniformOutput', false),'VerticalAlignment','bottom','HorizontalAlignment','right')
%         hold on, plot3(newTetCentroid(1), newTetCentroid(2), newTetCentroid(3), 'rx');
        isConvex = false;
        return
    end
end

% %Combine the two new tets created
% uniqueNodes = unique([nodeTets; Tnew]);
% %Obtain the CONCAVE shape connecting the tets (this is to allow concave
% %shapes between different tets)
% shapeBoundary = boundary(vertcat(Geo.Cells(uniqueNodes).X), 1);
% uniqueConnectedNodes = uniqueNodes(shapeBoundary);
% %Obtain how they are connected
% neighbours = unique(uniqueConnectedNodes(any(ismember(uniqueConnectedNodes, nodeItShouldBeConnected), 2), :));
% 
% 
% %Visualize if it is correct
% % allNodes = vertcat(Geo.Cells.X);
% % figure, trisurf(uniqueConnectedNodes, allNodes(:, 1), allNodes(:, 2), allNodes(:, 3))
% % text(allNodes(uniqueNodes, 1), allNodes(uniqueNodes, 2), allNodes(uniqueNodes, 3), cellfun(@num2str, num2cell(uniqueNodes), 'UniformOutput', false),'VerticalAlignment','bottom','HorizontalAlignment','right')
% 
% %It would be convex if the connections correspond to their assigned node
% %(nodeTets)
% isConvex = all(ismember(nodeTets, neighbours));
end

