function [Geo_0, Geo_n, Geo, Dofs, newYgIds, hasConverged] = PostFlip(Tnew, Ynew, oldTets, Geo, Geo_n, Geo_0, Dofs, newYgIds, Set, flipName, segmentToChange)
%POSTFLIP Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n; Geo_0_backup = Geo_0; Dofs_backup = Dofs;

Geo.log = sprintf('%s =>> %s-Flip: %i %i.\n', Geo.log, flipName, segmentToChange(1), segmentToChange(2));

[Geo] = RemoveTetrahedra(Geo, oldTets);
[Geo] = AddTetrahedra(Geo, Tnew, Ynew, Set);
[Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
[Geo_n] = AddTetrahedra(Geo_n, Tnew, Ynew, Set);
[Geo_0] = RemoveTetrahedra(Geo_0, oldTets);
[Geo_0] = AddTetrahedra(Geo_0, Tnew, Ynew, Set);

%try
    Geo   = Rebuild(Geo, Set);
    Geo   = BuildGlobalIds(Geo);
    Geo   = UpdateMeasures(Geo);
    
    Geo_n = Rebuild(Geo_n, Set);
    Geo_n = BuildGlobalIds(Geo_n);
    Geo_n = UpdateMeasures(Geo_n);
    
    Geo_0 = Rebuild(Geo_0, Set);
    Geo_0 = BuildGlobalIds(Geo_0);
        
    %% Update Geo_0 to be reset the vertices that we have changed averaging with previous Geo_0 and current Geo
    percentageGeo = 1 - Set.Reset_PercentageGeo0;
    for c=1:Geo.nCells
        if ismember(c, Tnew) && ~isempty(Geo.Cells(c).AliveStatus) && Geo.Cells(c).AliveStatus == 1
            Geo_0.Cells(c).X = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).X + percentageGeo * Geo.Cells(c).X;
            Geo_0.Cells(c).Y = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).Y + percentageGeo * Geo.Cells(c).Y;
            
            for f=1:length(Geo.Cells(c).Faces)
                Geo_0.Cells(c).Faces(f).Centre = Set.Reset_PercentageGeo0 * Geo_0.Cells(c).Faces(f).Centre + percentageGeo * Geo.Cells(c).Faces(f).Centre;
            end
        end
    end
    
    Geo_0 = UpdateMeasures(Geo_0);
    
% catch MException
%     Geo   = Geo_backup;
%     Geo_n = Geo_n_backup;
%     Geo_0 = Geo_0_backup;
%     Geo.log = strcat(Geo.log, sprintf('=>> %s-Flip rejected: ', flipName);
%     fprintf(MException.identifier);
%     fprintf('\n');
%     fprintf(MException.message);
%     fprintf('\n');
%     return
% end

if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
    PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
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
    
    newYgIds = unique([newYgIds; Geo.AssemblegIds]);
    Geo   = UpdateMeasures(Geo);
    Geo_n = UpdateMeasures(Geo_n);
    
    PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
    
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

