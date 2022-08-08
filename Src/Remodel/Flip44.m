function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip44(f, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n; Dofs_backup = Dofs;
Ys = Geo.Cells(numCell).Y;
Ts = Geo.Cells(numCell).T;
Face = Geo.Cells(numCell).Faces(f);

YsToChange=[Face.Tris(1).Edge(1); Face.Tris(2).Edge(1); Face.Tris(3).Edge(1); Face.Tris(4).Edge(1)];
[Ynew, Tnew] = YFlip44(Ys, Ts, YsToChange, Face, Geo);

targetTets = Geo.Cells(numCell).T(YsToChange,:);
Geo   = ReplaceYs(targetTets, Tnew, Ynew, Geo);
Geo_n = ReplaceYs(targetTets, Tnew, Ynew, Geo_n);

Geo   = RemoveFaces(f, Face.ij, Geo);
Geo_n = RemoveFaces(f, Face.ij, Geo_n);

%% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
Geo   = Rebuild(Geo, Set);
Geo_n = Rebuild(Geo_n, Set);

Geo   = BuildGlobalIds(Geo);
Geo_n = BuildGlobalIds(Geo_n);

Geo   = UpdateMeasures(Geo);
Geo_n = UpdateMeasures(Geo_n);

if ~CheckConvexity(Tnew,Geo_backup) && CheckTris(Geo)
    fprintf('=>> 44 Flip.\n');
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
        Dofs = Dofs_backup;
        fprintf('=>> 44-Flip rejected: did not converge\n');
        return
    end
    newYgIds = unique([newYgIds; Geo.AssemblegIds]);
    targetNodes = unique(targetTets);
    for n_i = 1:length(targetNodes)
        tNode = targetNodes(n_i);
        news = find(sum(ismember(Tnew,tNode)==1,2));
        if ~ismember(tNode, Geo.XgID)
            Geo_n.Cells(tNode).Y(end-length(news)+1:end,:) = Geo.Cells(tNode).Y(end-length(news)+1:end,:);
        end
    end
    Geo   = UpdateMeasures(Geo);
    Geo_n = UpdateMeasures(Geo_n);
    
    hasConverged = 1;
else
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    Dofs = Dofs_backup;
    fprintf('=>> 44-Flip rejected: is not compatible\n');
    return
end
end


