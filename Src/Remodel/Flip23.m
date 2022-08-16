function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip23(YsToChange, numCell, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n; Dofs_backup = Dofs;
Ys = Geo.Cells(numCell).Y;
Ts = Geo.Cells(numCell).T;

[Ynew, Tnew] = YFlip23(Ys, Ts, YsToChange, Geo);

ghostNodes = ismember(Tnew, Geo.XgID);
ghostNodes = all(ghostNodes, 2);
if any(ghostNodes)
    fprintf('=>> Flips 2-2 are not allowed for now\n');
    return
end

targetTets = Geo.Cells(numCell).T(YsToChange,:);

flipName = '2-3';
if isempty(Tnew) || CheckOverlappingTets(targetTets, Tnew, Geo, 'Internal')
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> %s-Flip rejected: is not compatible\n', flipName);
    return
end

Geo   = ReplaceYs(targetTets, Tnew, Ynew, Geo);
Geo_n = ReplaceYs(targetTets, Tnew, Ynew, Geo_n);

Geo   = Rebuild(Geo, Set, Tnew);
Geo_n = Rebuild(Geo_n, Set, Tnew);

Geo   = BuildGlobalIds(Geo);
Geo_n = BuildGlobalIds(Geo_n);

Geo   = UpdateMeasures(Geo);
Geo_n = UpdateMeasures(Geo_n);

if ~CheckConvexity(Tnew, Geo_backup) && CheckTris(Geo)
    fprintf('=>> 23 Flip.\n');
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
        Dofs = Dofs_backup;
        fprintf('=>> 23-Flip rejected: did not converge\n');
        return
    end
    newYgIds = unique([newYgIds; Geo.AssemblegIds]);
    Geo   = UpdateMeasures(Geo);
    Geo_n = UpdateMeasures(Geo_n);
    
    PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
    
    hasConverged = 1;
else
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    Dofs = Dofs_backup;
    fprintf('=>> 23-Flip rejected: is not compatible\n');
    return
end
end