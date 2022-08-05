function [Geo_n, Geo, Dofs, Set, newYgIds, hasConverged] = FlipAddNode(surroundingNodes, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIPADDNODE Summary of this function goes here
%   Detailed explanation goes here

hasConverged = 0;
Geo_backup = Geo; Geo_n_backup = Geo_n;

oldTets = Geo.Cells(surroundingNodes).T;
mainNode = surroundingNodes(~cellfun(@isempty, {Geo.Cells(surroundingNodes).AliveStatus}));
commonNodes = surroundingNodes;
commonNodes(ismember(commonNodes, mainNode)) = [];

if length(mainNode) > 1
    return
end

[Geo, newNodeIDs] = AddNewNode(Geo, mean(vertcat(Geo.Cells(commonNodes).X)));
[Geo_n] = AddNewNode(Geo_n, mean(vertcat(Geo.Cells(commonNodes).X)));

nodesToChange = [unique(oldTets); newNodeIDs];
[Tnew] = ConnectTetrahedra(Geo, nodesToChange, oldTets, mainNode);

%figure, tetramesh(Tnew, vertcat(Geo.Cells.X));
%figure, tetramesh(oldTets, vertcat(Geo.Cells.X));


if isempty(Tnew) || CheckOverlappingTets(oldTets, Tnew, Geo)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> AddNode-Flip rejected: is not compatible\n');
    return
end

fprintf('=>> AddNode-Flip.\n');

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
        fprintf('=>> AddNode-Flip rejected: did not converge\n');
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
    fprintf('=>> AddNode-Flip rejected: is not compatible\n');
    return
end

end

