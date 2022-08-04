function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipRemoveNode(nodeToRemove, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPREMOVENODE Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n;

oldTets = Geo.Cells(nodeToRemove).T;
nodesToChange = getNodeNeighbours(Geo, nodeToRemove);
mainNode = nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus}));

if length(mainNode) > 1
    return
end

if length(nodesToChange) > 4
    Tnew = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X)));
else
    Tnew = nodesToChange';
end

if length(nodesToChange) > 4 && CheckOverlappingTets(oldTets, Tnew, Geo)
    nodesToChange(~cellfun(@isempty, {Geo.Cells(nodesToChange).AliveStatus})) = [];
    [~,score] = pca(vertcat(Geo.Cells(nodesToChange).X));
    DT = delaunayTriangulation(score(:, 1:2));
    Tnew = horzcat(ones(size(DT.ConnectivityList, 1), 1) * mainNode, nodesToChange(DT.ConnectivityList));
end

if isempty(Tnew) || CheckOverlappingTets(oldTets, Tnew, Geo)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> RemoveNode-Flip rejected: is not compatible\n');
    return
end

fprintf('=>> RemoveNode-Flip.\n');

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
%% ----------------------------

if CheckTris(Geo) %%&& ~CheckConvexity(Tnew,Geo_backup)
    Dofs = GetDOFs(Geo, Set);
    [Dofs, Geo]  = GetRemodelDOFs(Tnew, Dofs, Geo);
    [Geo, Set, DidNotConverge] = SolveRemodelingStep(Geo_0, Geo_n, Geo, Dofs, Set);
    if DidNotConverge
        Geo   = Geo_backup;
        Geo_n = Geo_n_backup;
        fprintf('=>> RemoveNode-Flip rejected: did not converge\n');
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
    fprintf('=>> RemoveNode-Flip rejected: is not compatible\n');
    return
end
end

