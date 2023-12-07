function [surfacePoints] = getSurfacePointsCells(Geo, T)
%GETSURFACEPOINTSCELLS Summary of this function goes here
%   Detailed explanation goes here
precision = 0.1;
surfacePoints = {};
for tCell = T
    if ~isempty(Geo.Cells(tCell).AliveStatus) && Geo.Cells(tCell).AliveStatus == 1 && ~ismember(tCell, Geo.XgBottom)
        points = Geo.Cells(tCell).Y;
        regionIds = ones(size(points, 1), 1) == 1;
%         if any(ismember(T, Geo.XgTop))
%             regionIds = any(ismember(Geo.Cells(tCell).T, Geo.XgTop), 2);
%         elseif any(ismember(T, Geo.XgBottom))
%             regionIds = any(ismember(Geo.Cells(tCell).T, Geo.XgBottom), 2);
%         end
        x = points(regionIds, 1);
        y = points(regionIds, 2);
        z = points(regionIds, 3);

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
        surfacePoints(end+1) = {[x, y, z]};
    end
end
end

