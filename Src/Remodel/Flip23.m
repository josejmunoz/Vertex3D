function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)

hasConverged = 0;
Ys = Geo.Cells(numCell).Y;
Ts = Geo.Cells(numCell).T;
tetsToChange = Geo.Cells(numCell).T(YsToChange,:);
try
[Ynew, Tnew] = YFlip23(Ys, Ts, YsToChange, Geo);
catch
   return 
end

ghostNodes = ismember(Tnew, Geo.XgID);
ghostNodes = all(ghostNodes, 2);
if any(ghostNodes)
    fprintf('=>> Flips 2-2 are not allowed for now\n');
    return
end

% Rebuild topology and run mechanics
[Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, tetsToChange, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, 'Internal-23');
end