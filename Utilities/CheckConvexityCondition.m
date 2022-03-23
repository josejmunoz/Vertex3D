function [isConvex, tetID]=CheckConvexityCondition(Tnew,Tetrahedra, X)
%CHECKCONVEXITYCONDITION Summary of this function goes here
%   Check if the tetrahedron:
%   - is already created
%   - overlap with other tetrahedra
%   - is convex

isConvex = true;
tetID = -1;

%% Checking if the same tetrahadron is already on T
[foundTets, tetFoundIds] = ismember(sort(Tnew, 2),sort(Tetrahedra.DataRow, 2), 'rows');
if any(foundTets>0)
    tetID = tetFoundIds(foundTets);
    isConvex = false;
    return
end

%% Checking if Tnew overlap with other tetrahedra
% Here, we would calculate the plane of each face of the Tetrahedron and
% see if that intersect with other plane of a neighbouring tetrahedron. 
% Two faces would interesect if:
%  - They are not parallel AND
%  - They share a different point than the vertices on the corners
% %% Checking if Tnew overlap with other tetrahedra
% for numTnew = 1:size(Tnew, 1)
%     currentTet = Tnew(numTnew, :);
%     tetShape = alphaShape(X(currentTet, 1), X(currentTet, 2), X(currentTet, 3));
%     allXsExceptCurrentTet = 1:size(X, 1);
%     allXsExceptCurrentTet(Tnew(numTnew, :)) = [];
%     % Checking if any point of the Xs are inside the tetrahedra
%     if any(tetShape.inShape(X(allXsExceptCurrentTet, 1), X(allXsExceptCurrentTet, 2), X(allXsExceptCurrentTet, 3)))
%         tetID = numTnew;
%         isConvex = true;
%         return
%     end
% end

% for numTnew = 1:size(Tnew, 1)
%     currentTet = Tnew(numTnew, :);
%     tetShape = alphaShape(X(currentTet, 1), X(currentTet, 2), X(currentTet, 3));
%     allXsExceptCurrentTet = 1:size(X, 1);
%     allXsExceptCurrentTet(Tnew(numTnew, :)) = [];
%     % Checking if any point of the Xs are inside the tetrahedra
%     if any(tetShape.inShape(X(allXsExceptCurrentTet, 1), X(allXsExceptCurrentTet, 2), X(allXsExceptCurrentTet, 3)))
%         tetID = numTnew;
%         isConvex = true;
%         return
%     end
% end

%% Checking if Tnew is convex
disp('h')
for numTet = 1:size(Tnew, 1)
    currentTet = Tnew(numTet, :);
    for nextNumTet = numTet+1:size(Tnew, 1)
        nextTet = Tnew(nextNumTet, :);
        
        sharedVertices = currentTet(ismember(currentTet, nextTet));
        endPointCurrentTet = currentTet(ismember(currentTet, nextTet) == 0);
        endPointNextTet = nextTet(ismember(nextTet, currentTet) == 0);
        
        allPairedVertices = nchoosek(sharedVertices,2);
        for numPair = 1:size(allPairedVertices, 1)
            edge = allPairedVertices(numPair, :);
            
            theOtherVertex = sharedVertices(ismember(sharedVertices, edge) == 0);
            vectorToOtherVertex_Current = X(endPointCurrentTet, :) - X(theOtherVertex, :);
            vectorToOtherVertex_Next = X(endPointNextTet, :) - X(theOtherVertex, :);
            
%             figure, plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
%             hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
%             hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
            
            %https://stackoverflow.com/questions/2142552/calculate-the-angle-between-two-triangles-in-cuda
            % Here we create two normals of the triangles cross of two
            % vectors represents its normal
            normal2Triangle = cross(X(endPointNextTet, :) - X(edge(1), :), X(edge(2), :) - X(edge(1), :));
            normal1Triangle = cross(X(endPointCurrentTet, :) - X(edge(1), :), X(edge(2), :) - X(edge(1), :));
            
            if dot(vectorToOtherVertex_Current, normal1Triangle) / (norm(vectorToOtherVertex_Current) *  norm(normal1Triangle)) > 0 
                normal1Triangle = -normal1Triangle;
            end
            
            if dot(vectorToOtherVertex_Next, normal2Triangle) / (norm(vectorToOtherVertex_Next) *  norm(normal2Triangle)) > 0 
                normal2Triangle = -normal2Triangle;
            end
            
            % Define
            directionVectorTriangles = cross(normal1Triangle, normal2Triangle);
            directionVectorOrigin = cross(vectorToOtherVertex_Current, vectorToOtherVertex_Next);
            % Dot is how aligned are them: 1 is parallel; -1 opposite
            % direction
            alignmentDirections = dot(directionVectorOrigin, directionVectorTriangles) / (norm(directionVectorOrigin) * norm(directionVectorTriangles));
            
            determinant = det(directionVectorOrigin * directionVectorTriangles);
            atan2(determinant, alignmentDirections)
            
            %Another option is to do a dimensionality reduction and obtain
            %a 2D problem: https://math.stackexchange.com/questions/878785/how-to-find-an-angle-in-range0-360-between-2-vectors
            
            %
            if alignmentDirections < -0.5
                isConvex = false;
                return
            end
        end
    end
end

end

