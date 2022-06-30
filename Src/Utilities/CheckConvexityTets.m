function [isConvex] = CheckConvexityTets(nodeTets, Tnew, Geo)
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
end

