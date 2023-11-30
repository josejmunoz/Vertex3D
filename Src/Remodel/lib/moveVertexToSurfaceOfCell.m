function extrapolated_point_on_surface = moveVertexToSurfaceOfCell(Geo, T, newY)
    extrapolated_point_on_surface = [];
    precision = 0.1;

%     region
%     if any(ismember(T, Geo.XgTop))
%         
%     elseif any(ismember(T, Geo.XgBottom))
% 
%     end

    for tCell = T
        if ~isempty(Geo.Cells(tCell).AliveStatus) && Geo.Cells(tCell).AliveStatus == 1 && ~ismember(tCell, Geo.XgBottom)
            points = Geo.Cells(tCell).Y;
            x = points(:, 1);
            y = points(:, 2);
            z = points(:, 3);
            % Number of smoothing iterations
            shp = alphaShape(x, y, z, 1);
            pc = criticalAlpha(shp,'one-region');
            shp.Alpha = pc+3;
            [X, Y, Z] = meshgrid(min(x):precision:max(x), min(y):precision:max(y), min(z):precision:max(z));
            newVertices = [X(:), Y(:), Z(:)];
            inside_points = inShape(shp, X(:), Y(:), Z(:));
            newVertices = newVertices(inside_points, :);
            K = convhull(newVertices);
        
            x = newVertices(K, 1);
            y = newVertices(K, 2);
            z = newVertices(K, 3);
            
            % Point to be extrapolated
            extrapolation_point = newY;
            
            % Find the index of the closest point on the surface
            distances = sqrt((x - extrapolation_point(1)).^2 + (y - extrapolation_point(2)).^2 + (z - extrapolation_point(3)).^2);
            [~, min_index] = min(distances);
            
            % Extrapolated point on the surface
            extrapolated_point_on_surface(end+1, :) = [x(min_index), y(min_index), z(min_index)];
        end
    end
    extrapolated_point_on_surface = mean(extrapolated_point_on_surface);
end