function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip13(numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP13 Summary of this function goes here
%   Detailed explanation goes here
Geo_backup = Geo; Geo_n_backup = Geo_n; Dofs_backup = Dofs;
hasConverged = 0;

[~, numNeighbours_1, tetsNeighbours_1] = getVertexNeighbours(Geo, trisToChange.Edge(1), numCell);
[~, numNeighbours_2, tetsNeighbours_2] = getVertexNeighbours(Geo, trisToChange.Edge(2), numCell);

[trisArea_1, ~] = getTrisAreaOfNeighbours(Geo, Geo.Cells(numCell).globalIds(trisToChange.Edge(1)));
[trisArea_2, ~] = getTrisAreaOfNeighbours(Geo, Geo.Cells(numCell).globalIds(trisToChange.Edge(2)));

if ~all(trisArea_1 > Set.upperAreaThreshold)
    trisArea_1 = 0;
end

if ~all(trisArea_2 > Set.upperAreaThreshold)
    trisArea_2 = 0;
end
    
if isequal(trisArea_2, 0) && isequal(trisArea_1, 0)
    %fprintf('=>> 13-Flip rejected: not big enough triangles\n');
    return
end

if sum(trisArea_1 > Set.upperAreaThreshold) > sum(trisArea_2 > Set.upperAreaThreshold)
    vertexToExpand = trisToChange.Edge(1);
    oldTets = Geo.Cells(numCell).T(trisToChange.Edge(1), :);
else
    vertexToExpand = trisToChange.Edge(2);
    oldTets = Geo.Cells(numCell).T(trisToChange.Edge(2), :);
end

fprintf('=>> 13 Flip.\n');
nodesWithoutTheCell = Geo.Cells(numCell).T(vertexToExpand, :);
nodesWithoutTheCell(~cellfun(@isempty, {Geo.Cells(nodesWithoutTheCell).AliveStatus})) = [];
nodesWithoutCellPos = mean(vertcat(Geo.Cells(nodesWithoutTheCell).X));
newNodePosition = nodesWithoutCellPos;

[Geo, newNodeIDs] = AddNewNode(Geo, newNodePosition);
[Geo_n] = AddNewNode(Geo_n, newNodePosition);

[Tnew] = DoFlip13(Geo, newNodeIDs, oldTets');

if isempty(Tnew) || CheckOverlappingTets(oldTets, Tnew, Geo)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    Dofs = Dofs_backup;
    fprintf('=>> 13-Flip rejected: is not compatible\n');
    return
end

[Geo] = RemoveTetrahedra(Geo, oldTets);
[Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
[Geo] = AddTetrahedra(Geo, Tnew, Set);
[Geo_n] = AddTetrahedra(Geo_n, Tnew, Set);

%% All this, goes together when remodel occurs. TODO: PUT TOGETHER AS A FUNCTION
Geo   = Rebuild(Geo, Set);
Geo_n = Rebuild(Geo_n, Set);

Geo   = BuildGlobalIds(Geo);
Geo_n = BuildGlobalIds(Geo_n);

Geo   = UpdateMeasures(Geo);
Geo_n = UpdateMeasures(Geo_n);

if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup; 
        Dofs = Dofs_backup;
        fprintf('=>> 13-Flip rejected: did not converge\n');
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
    fprintf('=>> 13-Flip rejected: is not compatible\n');
    return
end
end

