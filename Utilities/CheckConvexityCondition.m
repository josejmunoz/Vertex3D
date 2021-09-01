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
        
        allPairedVertices = nchoosek(sharedVertices,2);        
        for numPair = 1:size(allPairedVertices, 1)
            edge = allPairedVertices(numPair, :);
            
            theOtherVertex = sharedVertices(ismember(sharedVertices, edge) == 0);

            if det([X(endPointCurrentTet, :); X(endPointNextTet, :); X(edge(2), :) - X(edge(1), :)]) < 0
                figure, plot3(X([endPointCurrentTet endPointNextTet], 1), X([endPointCurrentTet endPointNextTet], 2), X([endPointCurrentTet endPointNextTet], 3), 'x')
                hold on, plot3(X(edge, 1), X(edge, 2), X(edge, 3), 'bo')
                hold on, plot3(X(theOtherVertex, 1), X(theOtherVertex, 2), X(theOtherVertex, 3), 'ro')
                isConvex = false;
                tetID = [tetID; numTet, nextNumTet];
            end
        end
    end
end

if isConvex
    disp('All tetrahedra is convex');
elseif isequal(Tnew, Tetrahedra)
    %% Need to convexify
end

end

