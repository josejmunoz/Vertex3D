function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip24(Face, numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP24 Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n;

tetsToExpand = Geo.Cells(numCell).T(Face.Tris(trisToChange).Edge, :);
commonNodes = intersect(tetsToExpand(1, :), tetsToExpand(2, :));

fprintf('=>> 24 Flip.\n');
%% Pick the Ghost node
if ~isempty(firstNodeAlive)
    mainNode = Face.ij(1);
    commonNodes(commonNodes == Face.ij(1)) = [];
else
    mainNode = Face.ij(2);
    commonNodes(commonNodes == Face.ij(2)) = [];
end

[Geo, newNodeIDs] = AddNewNode(Geo, mean(vertcat(Geo.Cells(commonNodes).X)));
[Geo_n] = AddNewNode(Geo_n, mean(vertcat(Geo.Cells(commonNodes).X)));

%% Assign nodes to tets
oldTets = tetsToExpand;

[newTets] = ConnectTetrahedra(Geo, newNodeIDs, oldTets);

if isempty(Tnew) || CheckOverlappingTets(oldTets, Tnew, Geo)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> 24-Flip rejected: is not compatible\n');
    return
end

[Geo] = RemoveTetrahedra(Geo, oldTets);
[Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
[Geo] = AddTetrahedra(Geo, newTets, Set);
[Geo_n] = AddTetrahedra(Geo_n, newTets, Set);

%% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
Geo   = Rebuild(Geo, Set);
Geo_n = Rebuild(Geo_n, Set);

Geo   = BuildGlobalIds(Geo);
Geo_n = BuildGlobalIds(Geo_n);

Geo   = UpdateMeasures(Geo);
Geo_n = UpdateMeasures(Geo_n);
%% ----------------------------

%targetTets = testToSubstitute;
if CheckTris(Geo)
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(newTets, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
        fprintf('=>> 24-Flip rejected: did not converge\n');
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
    fprintf('=>> 24-Flip rejected: is not compatible\n');
    return
end
end

