function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip30(numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP30 Summary of this function goes here
%   Remove 1 node by combining two

    tetsToShrink = Geo.Cells(numCell).T(Face.Tris(trisToChange).Edge, :);
    commonNodes = intersect(tetsToShrink(1, :), tetsToShrink(2, :));
    opposingNodes = setxor(tetsToShrink(1, :), tetsToShrink(2, :));
    firstNodeAlive = Geo.Cells(Face.ij(1)).AliveStatus;
    secondNodeAlive = Geo.Cells(Face.ij(2)).AliveStatus;
    if xor(isempty(firstNodeAlive), isempty(secondNodeAlive))
        fprintf('=>> 30 Flip.\n');
        %% Pick the Ghost node
        if ~isempty(firstNodeAlive)
            commonNodes(commonNodes == Face.ij(1)) = [];
        else
            commonNodes(commonNodes == Face.ij(2)) = [];
        end
        
        % Check which tets overlap between the two 'commonNodes'
        [smallestDistance, commonNodeSmallest] = pdist2(vertcat(Geo.Cells(commonNodes).X), vertcat(Geo.Cells(opposingNodes).X), 'euclidean', 'Smallest', 1);
        [~, idsSmallestDistance]= min(smallestDistance);
        
        nodesToSubstitute = [commonNodes(commonNodeSmallest(idsSmallestDistance)), opposingNodes(idsSmallestDistance)];
        
        oldTets = vertcat(Geo.Cells(:).T);
        testToSubstitute = unique(sort(oldTets(sum(ismember(vertcat(Geo.Cells(:).T), nodesToSubstitute), 2) > 1, :), 2), 'row');
        
        [Geo, Tnew] = CombineTwoGhostNodes(Geo, Set, nodesToSubstitute);
        
        if isempty(Tnew)
            Geo   = Geo_backup;
            Geo_n = Geo_n_backup;
            fprintf('=>> 30-Flip rejected: is not compatible\n');
        end
        
        [Geo_n] = CombineTwoGhostNodes(Geo_n, Set, nodesToSubstitute);
    else
        return
    end
    
    %% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
    Geo   = Rebuild(Geo, Set);
    Geo_n = Rebuild(Geo_n, Set);
    
    Geo   = BuildGlobalIds(Geo);
    Geo_n = BuildGlobalIds(Geo_n);
    
    Geo   = UpdateMeasures(Geo);
    Geo_n = UpdateMeasures(Geo_n);
    %% ----------------------------
    
    if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
        Dofs = GetDOFs(Geo, Set);
        [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
        [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
        if DidNotConverge
            Geo   = Geo_backup;
            Geo_n = Geo_n_backup;
            fprintf('=>> 30-Flip rejected: did not converge\n');
        end
        
        newYgIds = unique([newYgIds; Geo.AssemblegIds]);
        Geo   = UpdateMeasures(Geo);
        Geo_n = UpdateMeasures(Geo_n);
        
        PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
        
        hasConverged = 1;
    else
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
        fprintf('=>> 30-Flip rejected: is not compatible\n');
    end
end

