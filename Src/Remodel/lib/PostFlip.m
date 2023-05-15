function [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName, segmentToChange)
%POSTFLIP Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n; Geo_0_backup = Geo_0; Dofs_backup = Dofs;

Geo.log = sprintf('%s =>> %s-Flip: %i %i.\n', Geo.log, flipName, segmentToChange(1), segmentToChange(2));

[Geo] = AddAndRebuildCells(Geo, oldTets, Tnew, Ynew, Set, 1);
Geo_n = Geo;
[Geo_0] = AddAndRebuildCells(Geo_0, oldTets, Tnew, Ynew, Set, 0);
%PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1)
%PostProcessingVTK(Geo_0, Geo_0, Set, Set.iIncr+2)

if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
    %PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+1)
    if Set.NeedToConverge
        Dofs = GetDOFs(Geo, Set);
        [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
        [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
        if DidNotConverge
            Geo   = Geo_backup;
            Geo_n = Geo_n_backup;
            Geo_0 = Geo_0_backup;
            Dofs = Dofs_backup;
            Geo.log = sprintf('%s =>> %s-Flip rejected: did not converge\n', Geo.log, flipName);
            return
        end
    end
    
    newYgIds = unique([newYgIds; Geo.AssemblegIds]);
    Geo   = UpdateMeasures(Geo);
    
    hasConverged = 1;
else
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    Dofs = Dofs_backup;
    Geo_0 = Geo_0_backup;
    Geo.log = sprintf('%s =>> %s-Flip rejected: is not compatible\n', Geo.log, flipName);
    return
end
end

