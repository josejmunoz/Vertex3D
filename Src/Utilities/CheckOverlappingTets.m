function [overlaps, correctTets] = CheckOverlappingTets(oldTets, newTets, Geo, flipType)
%CHECKOVERLAPPINGTETS Summary of this function goes here
%   Detailed explanation goes here

overlaps = 0;
correctTets = [];

newVol = 0;
volumes = [];
for tet = newTets'
    [vol] = ComputeTetVolume(tet, Geo);
    volumes(end+1) = vol;
    newVol = newVol + vol;
end

normVols = volumes/max(volumes);
if any(normVols < 0.07)
    correctTets = newTets(normVols >= 0.07, :);
    overlaps = 1;
    return;
end


if isequal(flipType, 'AddNode')
    
elseif isequal(flipType, 'RemoveNode')
    %% Check nodes connectivity
    nodesToCheck = intersect(unique(oldTets), unique(newTets));
    nodesToCheck = nodesToCheck(cellfun(@isempty, {Geo.Cells(nodesToCheck).AliveStatus}));

    [Geo] = RemoveTetrahedra(Geo, oldTets);
    [Geo] = AddTetrahedra(Geo, newTets);
    newNeighs = arrayfun(@(x) getNodeNeighbours(Geo, x), nodesToCheck, 'UniformOutput', false);
    newNeighs = cellfun(@(x) x(ismember(x, nodesToCheck)), newNeighs, 'UniformOutput', false);

    if sum(cellfun(@length, newNeighs))/2 > (2*length(nodesToCheck) - 3)
        overlaps = 1;
        return
    end
elseif isequal(flipType, 'Internal')
    %% Check if the volume from previous space is the same occupied by the new tets
    oldVol = 0;
    for tet = oldTets'
        [vol] = ComputeTetVolume(tet, Geo);
        oldVol = oldVol + vol;
    end
    
    if abs(newVol - oldVol)/newVol > 0.05
        overlaps = 1;
        return
    end
%     %% Check if some skinny tet has been created
%     if any(CheckSkinnyTets(newTets, Geo))
%         overlaps = 1;
%         return
%     end
end

%% Check if they overlap
for numTet = 1:size(newTets, 1)-1
    currentTet = newTets(numTet, :);
    for nextNumTet = numTet+1:size(newTets, 1)
        nextTet = newTets(nextNumTet, :);
        if ~isequal(currentTet, nextTet)
            %% 1st Shape
            shape1 = vertcat(Geo.Cells(currentTet).X);
            
            %% 2nd Shape
            shape2 = vertcat(Geo.Cells(nextTet).X);
            reorderingTet = delaunayTriangulation(shape2);
            shape2 = shape2(reorderingTet.ConnectivityList, :);
            
            %             %% GJK mws262: Check if they overlap
            %             % There are some faults with winter.dev conversion hence using GJK
            %             % Collision Detection from mws262
            %             % https://github.com/mws262/MATLAB-GJK-Collision-Detection
            %             %Point 1 and 2 selection (line segment)
            %             direction = [1 0 0];
            %             [points] = simplex_line(direction, shape2, shape1);
            %
            %             %Point 3 selection (triangle)
            %             [points,overlaps] = simplex_triangle(points, shape2, shape1);
            %
            %             %Point 4 selection (tetrahedron)
            %             if overlaps == 1 %Only bother if we could find a viable triangle.
            %                 [points, overlaps] = simplex_tetrahedron(points, shape2, shape1);
            %             end
            
            %% https://uk.mathworks.com/matlabcentral/answers/327990-generate-random-coordinates-inside-a-convex-polytope
            %tic
            shape1 = vertcat(Geo.Cells(currentTet).X);
            CH = convhull(shape1);
            ntri = size(CH, 1);
            xycent = mean(shape1,1);
            nxy = size(shape1,1);
            ncent = nxy+1;
            shape1(ncent,:) = xycent;
            tri = [CH,repmat(ncent,ntri,1)];
            
            xy = shape1;
            %             figure
            %             plot3(xy(:,1),xy(:,2), xy(:,3),'bo');
            %             hold on
            %             plot3([xy(tri(:,1),1),xy(tri(:,2),1),xy(tri(:,3),1), xy(tri(:,4),1)]',[xy(tri(:,1),2),xy(tri(:,2),2),xy(tri(:,3),2),xy(tri(:,4),2)]',[xy(tri(:,1),3),xy(tri(:,2),3),xy(tri(:,3),3),xy(tri(:,4),3)]','g-')
            
            V = zeros(1,ntri);
            for ii = 1:ntri
                V(ii) = abs(det(xy(tri(ii,1:3),:) - xycent));
            end
            V = V/sum(V);
            M = 1000;
            [~,~,simpind] = histcounts(rand(M,1),cumsum([0,V]));
            
            r1 = rand(M,1);
            uvw = xy(tri(simpind,1),:).*r1 + xy(tri(simpind,2),:).*(1-r1);
            r2 = sqrt(rand(M,1));
            uvw = uvw.*r2 + xy(tri(simpind,3),:).*(1-r2);
            r3 = nthroot(rand(M,1),3);
            uvw = uvw.*r3 + xy(tri(simpind,4),:).*(1-r3);
            %plot3(uvw(:,1),uvw(:,2),uvw(:,3),'.')
            
            
            aShape2 = alphaShape(shape2);
            for numPoint = size(uvw, 1)
                if aShape2.inShape(uvw(numPoint, :))
                    overlaps = 1;
                    %tetramesh(vertcat(currentTet, nextTet), vertcat(Geo.Cells.X))
                    return
                end
            end
            %tetramesh(vertcat(currentTet, nextTet), vertcat(Geo.Cells.X))
            %toc
        end
    end
end

end

