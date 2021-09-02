function [isConvex, tetID]=CheckConvexityCondition(Tnew,Tetrahedra, X)
%CHECKCONVEXITYCONDITION Summary of this function goes here
%   Check if the tetrahedron:
%   - is already created
%   - overlap with other tetrahedra
%   - is convex

isConvex = true;
tetID = [];

%% Checking if the same tetrahadron is already on T
if isempty(Tnew) == 0
    [foundTets, tetFoundIds] = ismember(sort(Tnew, 2),sort(Tetrahedra.DataRow, 2), 'rows');
    if any(foundTets>0)
        tetID = tetFoundIds(foundTets);
        isConvex = false;
        return
    end
else
    Tnew = Tetrahedra.DataRow;
end

%% Checking if Tnew overlap with other tetrahedra
% Here, we would calculate the plane of each face of the Tetrahedron and
% see if that intersect with other plane of a neighbouring tetrahedron.
%X = X * 100;
for numTet = 1:size(Tnew, 1)
    currentTet = Tnew(numTet, :);
    for nextNumTet = numTet+1:size(Tnew, 1)
        nextTet = Tnew(nextNumTet, :);
        
        sharedVertices = currentTet(ismember(currentTet, nextTet));
        endPointCurrentTet = currentTet(ismember(currentTet, nextTet) == 0);
        endPointNextTet = nextTet(ismember(nextTet, currentTet) == 0);
        
        if length(sharedVertices) ~= 3 % They don't share a face
            continue
        end
        
        %In order to maintain the vertices order, we do it manually
        allPairedVertices = [sharedVertices(1), sharedVertices(2);
            sharedVertices(2), sharedVertices(3);
            sharedVertices(3), sharedVertices(1)];
        for numPair = 1:size(allPairedVertices, 1)
            edge = allPairedVertices(numPair, :);
            
            theOtherVertex = sharedVertices(ismember(sharedVertices, edge) == 0);
            
            tetVertices = [endPointCurrentTet endPointNextTet edge theOtherVertex];
            K = convhull(X(tetVertices, 1), X(tetVertices, 2), X(tetVertices, 3));
            K_sorted = sort(K, 2);
            if ismember([1 2], K_sorted(:, 1:2), 'rows') && all(sum(ismember(K, [3 4]), 2) < 2)
%                 h = figure;
%                 plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
%                 %hold on, quiver3(X(endPointNextTet, 1), X(endPointNextTet, 2), X(endPointNextTet, 3), normalNextTriangle(1), normalNextTriangle(2), normalNextTriangle(3));
%                 %hold on, quiver3(X(endPointCurrentTet, 1), X(endPointCurrentTet, 2), X(endPointCurrentTet, 3), normalCurrentTriangle(1), normalCurrentTriangle(2), normalCurrentTriangle(3));
%                 hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
%                 hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
%                 trisurf(K,X(tetVertices,1),X(tetVertices,2),X(tetVertices,3));
%                 a = alphaShape(X(tetVertices, 1), X(tetVertices, 2), X(tetVertices, 3));
%                 plot(a, 'FaceAlpha', 0.5)
%                 close(h)
                isConvex = false;
                tetID = [tetID; numTet, nextNumTet];
            end
            
            
%             vectorToOtherVertex_Current = X(endPointCurrentTet, :) - X(theOtherVertex, :);
%             vectorToOtherVertex_Next = X(endPointNextTet, :) - X(theOtherVertex, :);
% 
%             % Here we create two normals of the triangles cross of two
%             % vectors represents its normal
%             normalNextTriangle = cross(X(edge(2), :) - X(endPointNextTet, :), X(edge(1), :) - X(endPointNextTet, :));
%             normalCurrentTriangle = cross(X(edge(2), :) - X(endPointCurrentTet, :), X(edge(1), :) - X(endPointCurrentTet, :));
%             
%             if dot(vectorToOtherVertex_Current, normalCurrentTriangle) / (norm(vectorToOtherVertex_Current) *  norm(normalCurrentTriangle)) > 0 
%                 normalCurrentTriangle = -normalCurrentTriangle;
%             end
%             
%             if dot(vectorToOtherVertex_Next, normalNextTriangle) / (norm(vectorToOtherVertex_Next) *  norm(normalNextTriangle)) > 0 
%                 normalNextTriangle = -normalNextTriangle;
%             end
%             
%             h = figure;
%             plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
%             hold on, quiver3(X(endPointNextTet, 1), X(endPointNextTet, 2), X(endPointNextTet, 3), normalNextTriangle(1), normalNextTriangle(2), normalNextTriangle(3));
%             hold on, quiver3(X(endPointCurrentTet, 1), X(endPointCurrentTet, 2), X(endPointCurrentTet, 3), normalCurrentTriangle(1), normalCurrentTriangle(2), normalCurrentTriangle(3));
%             hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
%             hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
%             
%             parallelepidVectors = [normalCurrentTriangle; normalNextTriangle; X(edge(1), :) - X(edge(2), :)];
% %             minValueAllowed = 0.000001;
% %             minValues = min(parallelepidVectors);
% %             parallelepidVectors(:, minValues < minValueAllowed) = parallelepidVectors(:, minValues < minValueAllowed) + (minValueAllowed - minValues(minValues < minValueAllowed));
%             det(parallelepidVectors)
%             edge
%             close(h)
%             if det(parallelepidVectors) < 0
%                 h = figure;
%                 plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
%                 hold on, quiver3(X(endPointNextTet, 1), X(endPointNextTet, 2), X(endPointNextTet, 3), normalNextTriangle(1), normalNextTriangle(2), normalNextTriangle(3));
%                 hold on, quiver3(X(endPointCurrentTet, 1), X(endPointCurrentTet, 2), X(endPointCurrentTet, 3), normalCurrentTriangle(1), normalCurrentTriangle(2), normalCurrentTriangle(3));
%                 hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
%                 hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
%                 a = alphaShape(X([endPointCurrentTet endPointNextTet edge theOtherVertex], 1), X([endPointCurrentTet endPointNextTet edge theOtherVertex], 2), X([endPointCurrentTet endPointNextTet edge theOtherVertex], 3));
%                 plot(a)
%                 close(h)
%                 isConvex = false;
%                 tetID = [tetID; numTet, nextNumTet];
%             end
        end
    end
end

if isConvex
    disp('All tetrahedra is convex');
elseif isequal(Tnew, Tetrahedra)
    %% Need to convexify
end

end

