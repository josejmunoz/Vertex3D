function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip42(Face, numCell, trisToChange, Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP42 Summary of this function goes here
%   Remove 1 node

Geo_backup = Geo; Geo_n_backup = Geo_n;

tetsToShrink = Geo.Cells(numCell).T(Face.Tris(trisToChange).Edge, :);
commonNodes = intersect(tetsToShrink(1, :), tetsToShrink(2, :));

fprintf('=>> 42 Flip.\n');
%% Pick the Ghost node
if ~isempty(firstNodeAlive)
    mainNode = Face.ij(1);
    commonNodes(commonNodes == Face.ij(1)) = [];
else
    mainNode = Face.ij(2);
    commonNodes(commonNodes == Face.ij(2)) = [];
end

%Check commonNodes neighbourhood
commonNeighboursOfNodes = {getNodeNeighbours(Geo, commonNodes(1)); getNodeNeighbours(Geo, commonNodes(2))};

% Get the smallest neighbourhood between the two nodes
[~, smallestNumNeighs] = min([length(commonNeighboursOfNodes{1}), length(commonNeighboursOfNodes{2})]);

neighboursToUse = commonNeighboursOfNodes{smallestNumNeighs};
nodeToRemove = commonNodes(smallestNumNeighs);

neighboursToUse(neighboursToUse==mainNode) = [];

%Remove the selected node from that neighbourhood and
%reconnect them
oldTets = Geo.Cells(nodeToRemove).T;

nodeNeighbours = arrayfun(@(x) getNodeNeighbours(Geo, x), neighboursToUse, 'UniformOutput', false);
nodesConnected = neighboursToUse;
missingNeighbours = {};
for numNode = neighboursToUse'
    numNodePosition = find(neighboursToUse==numNode);
    missingNeighbours{numNodePosition} = neighboursToUse(~ismember(neighboursToUse, nodeNeighbours{neighboursToUse==numNode}));
    
    if ismember(commonNodes(setdiff(1:2, smallestNumNeighs)), missingNeighbours{numNodePosition})
        missingNeighbours{numNodePosition} = [];
    else
        nodesConnected = intersect(nodesConnected, missingNeighbours{numNodePosition});
    end
end

Tnew = [];
if length(nodesConnected) == 2
    opposingNodes = setdiff(neighboursToUse, nodesConnected);
    nodesToChange = [mainNode; nodesConnected; opposingNodes];
    Tnew = nodesToChange(delaunayn(vertcat(Geo.Cells(nodesToChange).X)));
else
    fprintf('NEED TO CHECKKKKK!!');
end

if isempty(Tnew)
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> 42-Flip rejected: is not compatible\n');
end

[Geo] = RemoveTetrahedra(Geo, oldTets);
[Geo_n] = RemoveTetrahedra(Geo_n, oldTets);
[Geo] = AddTetrahedra(Geo, Tnew, Set);
[Geo_n] = AddTetrahedra(Geo_n, Tnew, Set);


%% TODO: CHECK THIS!
%[overlaps] = CheckOverlappingTets(goodTets, testTets, Geo);

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
        fprintf('=>> 42-Flip rejected: did not converge\n');
    end
    
    newYgIds = unique([newYgIds; Geo.AssemblegIds]);
    Geo   = UpdateMeasures(Geo);
    Geo_n = UpdateMeasures(Geo_n);
    
    PostProcessingVTK(Geo, Geo_0, Set, Set.iIncr+2)
else
    Geo   = Geo_backup;
    Geo_n = Geo_n_backup;
    fprintf('=>> 42-Flip rejected: is not compatible\n');
end
end

