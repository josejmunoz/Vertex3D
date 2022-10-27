function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip44(f, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)

Ys = Geo.Cells(numCell).Y;
Ts = Geo.Cells(numCell).T;
Face = Geo.Cells(numCell).Faces(f);

YsToChange=[Face.Tris(1).Edge(1); Face.Tris(2).Edge(1); Face.Tris(3).Edge(1); Face.Tris(4).Edge(1)];
[~, Tnew] = YFlip44(Ys, Ts, YsToChange, Face, Geo);

tetsToChange = Geo.Cells(numCell).T(YsToChange,:);

% Rebuild topology and run mechanics
[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, tetsToChange, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'Internal-44');
end


