function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = Flip13(numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP13 Summary of this function goes here
%   Detailed explanation goes here
Geo_backup = Geo; Geo_n_backup = Geo_n;
hasConverged = 0;

[~, numNeighbours_1, tetsNeighbours_1] = getVertexNeighbours(Geo, trisToChange.Edge(1), numCell);
[~, numNeighbours_2, tetsNeighbours_2] = getVertexNeighbours(Geo, trisToChange.Edge(2), numCell);

if numNeighbours_1 == numNeighbours_2
    [trisArea_1, ~] = getTrisAreaOfNeighbours(Geo, Geo.Cells(numCell).globalIds(trisToChange.Edge(1)));
    [trisArea_2, ~] = getTrisAreaOfNeighbours(Geo, Geo.Cells(numCell).globalIds(trisToChange.Edge(2)));
    
    if sum(trisArea_1) > sum(trisArea_2)
        vertexToExpand = trisToChange.Edge(1);
        oldTets = Geo.Cells(numCell).T(trisToChange.Edge(1), :);
    else
        vertexToExpand = trisToChange.Edge(2);
        oldTets = Geo.Cells(numCell).T(trisToChange.Edge(2), :);
    end
elseif numNeighbours_1 < numNeighbours_2
    vertexToExpand = trisToChange.Edge(1);
    oldTets = Geo.Cells(numCell).T(trisToChange.Edge(1), :);
elseif numNeighbours_1 > numNeighbours_2
    vertexToExpand = trisToChange.Edge(2);
    oldTets = Geo.Cells(numCell).T(trisToChange.Edge(2), :);
end

fprintf('=>> 13 Flip.\n');
nodesWithoutTheCell = Geo.Cells(numCell).T(vertexToExpand, :);
nodesWithoutTheCell(~cellfun(@isempty, {Geo.Cells(nodesWithoutTheCell).AliveStatus})) = [];
allNodesPos = mean(vertcat(Geo.Cells(Geo.Cells(numCell).T(vertexToExpand, :)).X));
nodesWithoutCellPos = mean(vertcat(Geo.Cells(nodesWithoutTheCell).X));
newNodePosition = [allNodesPos(1:2), nodesWithoutCellPos(3)];

[Geo, newNodeIDs] = AddNewNode(Geo, newNodePosition);
[Geo_n] = AddNewNode(Geo_n, newNodePosition);

[newTets] = DoFlip13(Geo, newNodeIDs, oldTets');

if isempty(newTets)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> 13-Flip rejected: is not compatible\n');
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

if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(newTets, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
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
    fprintf('=>> 13-Flip rejected: is not compatible\n');
    return
end
end

