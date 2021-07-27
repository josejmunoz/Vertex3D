function [isConvex, tetID]=CheckConvexityCondition(Tnew,Tetrahedra, X)
%CHECKCONVEXITYCONDITION Summary of this function goes here
%   Check if the tetrahedron:
%   - is already created
%   - overlap with other tetrahedra
%   - is convex

isConvex = false;
tetID = -1;

%% Checking if the same tetrahadron is already on T
[foundTets, tetFoundIds] = ismember(sort(Tnew, 2),sort(Tetrahedra.DataRow, 2), 'rows');
if any(foundTets>0)
    tetID = tetFoundIds(foundTets);
    isConvex = true;
    return
end

%% Checking if Tnew overlap with other tetrahedra
% Here, we would calculate the plane of each face of the Tetrahedron and
% see if that intersect with other plane of a neighbouring tetrahedron. 
% Two faces would interesect if:
%  - They are not parallel AND
%  - They share a different point than the vertices on the corners

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

end

