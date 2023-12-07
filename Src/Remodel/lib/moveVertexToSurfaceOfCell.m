function Geo = moveVertexToSurfaceOfCell(Geo, T, surfacePoints)
    x = surfacePoints(:, 1);
    y = surfacePoints(:, 2);
    z = surfacePoints(:, 3);

    for tCell = T'
        if ~isempty(Geo.Cells(tCell).AliveStatus) && Geo.Cells(tCell).AliveStatus == 1 && ~ismember(tCell, Geo.XgBottom)
            newY = Geo.Cells(tCell).Y(ismember(sort(Geo.Cells(tCell).T, 2), sort(T)', 'rows'), :);
            % Find the index of the closest point on the surface
            distances = sqrt((x - newY(1)).^2 + (y - newY(2)).^2 + (z - newY(3)).^2);
            [~, min_index] = min(distances);
            
            % Extrapolated point on the surface
            Geo.Cells(tCell).Y(ismember(sort(Geo.Cells(tCell).T, 2), sort(T)', 'rows'), :) = [x(min_index), y(min_index), z(min_index)];
        end
    end
end