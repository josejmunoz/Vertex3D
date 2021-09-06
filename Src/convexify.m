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
        newXs(currentVertices(1:2), :) = newXs(currentVertices(1:2), :) - difference;
        
        Tnew = tetrahedra(sum(ismember(tetrahedra, currentVertices(1:2)), 2)>3, :);
        %Check if the current Tet is now convex
        if CheckConvexityCondition(Tnew, tetrahedra, newXs, false)
            % We go back to the old vertices position
            newXs(currentVertices(1:2), :) = X(currentVertices(1:2), :);
        end
        
        Tnew = tetrahedra(sum(ismember(tetrahedra, currentVertices(1:2)), 2)>0, :);
        [isConvex_new, concaveTetVertices_new] = CheckConvexityCondition(Tnew, tetrahedra, newXs, false);
        [isConvex_old, concaveTetVertices_old] = CheckConvexityCondition(Tnew, tetrahedra, X, false);
        if size(concaveTetVertices_new, 1) <= size(concaveTetVertices_old, 1)
            newXs(currentVertices(1:2), :) = X(currentVertices(1:2), :);
        else % Current tetrahedron is now convex
            continue;
        end
        %% Method 2: Move the other vertices 
        sum(ismember(K.ConnectivityList, tetrahedra, 'rows'))
        
        
    end
    disp('Convexify - done');
end

