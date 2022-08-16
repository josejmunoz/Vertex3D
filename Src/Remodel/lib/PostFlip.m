function [Geo, Geo_n, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName)
%POSTFLIP Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n; Dofs_backup = Dofs;

if isempty(Tnew) || CheckOverlappingTets(oldTets, Tnew, Geo, flipName)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> %s-Flip rejected: is not compatible\n', flipName);
    return
end

fprintf('=>> %s-Flip.\n', flipName);

[Geo] = RemoveTetrahedra(Geo, oldTets);
[Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
[Geo] = AddTetrahedra(Geo, Tnew, Set);
[Geo_n] = AddTetrahedra(Geo_n, Tnew, Set);

try
    Geo   = Rebuild(Geo, Set, Tnew);
    Geo_n = Rebuild(Geo_n, Set, Tnew);
    
    Geo   = BuildGlobalIds(Geo);
    Geo_n = BuildGlobalIds(Geo_n);
    
    Geo   = UpdateMeasures(Geo);
    Geo_n = UpdateMeasures(Geo_n);
catch MException
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> %s-Flip rejected: ', flipName);
    fprintf(MException.identifier);
    fprintf('\n');
    fprintf(MException.message);
    fprintf('\n');
    return
end

if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
        Dofs = Dofs_backup;
        fprintf('=>> %s-Flip rejected: did not converge\n', flipName);
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
    fprintf('=>> %s-Flip rejected: is not compatible\n', flipName);
    return
end
end

