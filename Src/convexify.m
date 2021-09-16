function [newXs] = convexify(tetrahedra, X, concaveTetVertices)
%CONVEXIFY Convexify tetrahedra
%   Detailed explanation goes here
    newXs = X;
    for numPairOfTets = 1:size(concaveTetVertices, 1)
        %Order is: [endPointCurrentTet endPointNextTet edgeVertices theOtherVertex]
        currentVertices = concaveTetVertices(numPairOfTets, :);
        
        %% Method 1: Increase the height of the 'edge vertices'
        outerVertices = mean(newXs(currentVertices(1:2), :));
        edgeVertices = mean(newXs(currentVertices(3:4), :));
        
        difference = outerVertices - edgeVertices;
        
        %Modify the outer vertices to make them a 
        indicesToCheck = 1:2;
        newXs(currentVertices(indicesToCheck), :) = newXs(currentVertices(indicesToCheck), :) - (difference*1.0001);
        
        [newXs, isConvex] = checkCurrentTetIsNowConvex(newXs, tetrahedra, currentVertices, X, indicesToCheck);
        
        if isConvex
            continue 
        end
        
        indicesToCheck = 3:4;
        newXs(currentVertices(indicesToCheck), :) = newXs(currentVertices(indicesToCheck), :) + (difference*1.0001);
        
        [newXs, isConvex] = checkCurrentTetIsNowConvex(newXs, tetrahedra, currentVertices, X, indicesToCheck);
        if isConvex
            continue 
        end
        %% Method 2: Move the other vertices 
        %sum(ismember(K.ConnectivityList, tetrahedra, 'rows'))
        numPairOfTets
    end
    [currentFace_isConvex, concaveTetVertices_new] = CheckConvexityCondition([], tetrahedra, newXs, false);
    disp('Convexify - done');
end

function [newXs, isConvex] = checkCurrentTetIsNowConvex(newXs, tetrahedra, currentVertices, X, indicesToCheck)
%% Check if the current Tet is now convex and has not added more concavities

isConvex = false;

Tnew = tetrahedra(sum(ismember(tetrahedra, currentVertices), 2)>3, :);
[currentFace_isConvex, concaveTetVertices] = CheckConvexityCondition(Tnew, tetrahedra, newXs, false);
if currentFace_isConvex == false
    % We go back to the old vertices position
    newXs(currentVertices(indicesToCheck), :) = X(currentVertices(indicesToCheck), :);
    return
end

Tnew = tetrahedra(sum(ismember(tetrahedra, currentVertices(indicesToCheck)), 2)>0, :);
[isConvex_new, concaveTetVertices_new] = CheckConvexityCondition(Tnew, tetrahedra, newXs, false);
[isConvex_old, concaveTetVertices_old] = CheckConvexityCondition(Tnew, tetrahedra, X, false);
if size(concaveTetVertices_new, 1) >= size(concaveTetVertices_old, 1)
    newXs(currentVertices(indicesToCheck), :) = X(currentVertices(indicesToCheck), :);
else % Current tetrahedron is now convex
    isConvex = true;
end

end

