function [overlaps] = CheckOverlappingTets(goodTets, testTets, Geo)
%CHECKOVERLAPPINGTETS Summary of this function goes here
%   Detailed explanation goes here

for numTet = 1:size(testTets, 1)-1
    currentTet = testTets(numTet, :);
    for nextNumTet = numTet+1:size(testTets, 1)
        nextTet = testTets(nextNumTet, :);
        if ~isequal(currentTet, nextTet)
            %% 1st Shape
            testXs = vertcat(Geo.Cells(currentTet).X);
            [faces] = convhull(testXs);
            shape1 = patch('Faces', faces, 'Vertices', testXs, 'visible','off');

            %% 2nd Shape
            testXs = vertcat(Geo.Cells(nextTet).X);
            [faces] = convhull(testXs);
            shape2 = patch('Faces', faces, 'Vertices', testXs, 'visible','off');
            
            %% Check if they overlap
            overlaps = GJK(shape1, shape2, 2);
            if overlaps
                return
            end
        end
    end
end

% 
% overlaps = true;
% for goodTet = 1:size(goodTets, 1)
%     for numTnew = 1:size(testTets, 1)
%         %Get the centre of the tetrahedron test tets and look for if it is inside any
%         %other tet.
%         newTet = testTets(numTnew, :);
%         newTetCentroid = mean(vertcat(Geo.Cells(newTet).X), 1);
% 
%         tetXs = vertcat(Geo.Cells(goodTets(goodTet)).X);
%         tetShape = alphaShape(tetXs);
% 
%         % Checking if any point of the Xs are inside the tetrahedra
%         if any(tetShape.inShape(newTetCentroid(1), newTetCentroid(2), newTetCentroid(3)))
%     %         figure, plot(tetShape)
%     %         text(tetXs(:, 1), tetXs(:, 2), tetXs(:, 3), cellfun(@num2str, num2cell(nodeTets), 'UniformOutput', false),'VerticalAlignment','bottom','HorizontalAlignment','right')
%     %         hold on, plot3(newTetCentroid(1), newTetCentroid(2), newTetCentroid(3), 'rx');
%             overlaps = false;
%             return
%         end
%     end
% end
end

