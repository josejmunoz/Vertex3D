function [Geo, Set] = InitializeGeometry_3DVoronoi(Geo, Set)
%INITIALIZEGEOMETRY_3DVORONOI Summary of this function goes here
%   Detailed explanation goes here

nCells = 50;
nSeeds = nCells*3;
lloydIterations = 100;
distorsion = 0;
cellHeight = 1;

rng default
x = rand(nSeeds, 1);
y = rand(nSeeds, 1);

seedsXY = horzcat(x,y);

%% Get central
distanceSeeds = pdist2(seedsXY, [0.5 0.5]);
[~, indices] = sort(distanceSeeds);
seedsXY = seedsXY(indices, :);

%% Homogeneize voronoi diagram
for numIter = 1:lloydIterations
    DT = delaunayTriangulation(seedsXY);
    [V, D] = voronoiDiagram(DT);
    for numCell = 1:nCells
        currentVertices = V(D{numCell}, :);
        seedsXY(numCell, :) = mean(currentVertices);
    end
end

%% Build 3D topology
[trianglesConnectivity, neighboursNetwork, cellEdges, verticesOfCell_pos] = Build3DVoronoiTopo(seedsXY);

seedsXY_topoChanged = [seedsXY(:, 1), seedsXY(:, 2) + distorsion];
% seedsXY_topoChanged(:, 2) = seedsXY_topoChanged(:, 2) - min(seedsXY_topoChanged(:, 2));
% seedsXY_topoChanged(:, 2) = seedsXY_topoChanged(:, 2) / max(seedsXY_topoChanged(:, 2));

[trianglesConnectivity_topoChanged, neighboursNetwork_topoChanged, cellEdges_topoChanged, verticesOfCell_pos_topoChanged] = Build3DVoronoiTopo(seedsXY_topoChanged);

%% Create node connections:
X(:, 1) = mean([seedsXY(:, 1), seedsXY_topoChanged(:, 1)], 2);
X(:, 2) = mean([seedsXY(:, 2), seedsXY_topoChanged(:, 2)], 2);
X(:, 3) = zeros(1, size(X, 1));
XgTopFaceCentre = horzcat(seedsXY, repmat(cellHeight, length(seedsXY), 1));
XgBottomFaceCentre = horzcat(seedsXY_topoChanged, repmat(-cellHeight, length(seedsXY_topoChanged), 1));
XgTopVertices = [verticesOfCell_pos, repmat(cellHeight, size(verticesOfCell_pos, 1), 1)];
XgBottomVertices = [verticesOfCell_pos_topoChanged, repmat(-cellHeight, size(verticesOfCell_pos_topoChanged, 1), 1)];

X_bottomNodes = vertcat(XgBottomFaceCentre, XgBottomVertices);
X_bottomIds = size(X, 1) + 1: size(X, 1) + size(X_bottomNodes, 1);
X_bottomFaceIds = X_bottomIds(1:size(XgBottomFaceCentre, 1));
X_bottomVerticesIds = X_bottomIds(size(XgBottomFaceCentre, 1)+1:end);
X = vertcat(X, X_bottomNodes);

X_topNodes = vertcat(XgTopFaceCentre, XgTopVertices);
X_topIds = size(X, 1) + 1: size(X, 1) + size(X_topNodes, 1);
X_topFaceIds = X_topIds(1:size(XgTopFaceCentre, 1));
X_topVerticesIds = X_topIds(size(XgTopFaceCentre, 1)+1:end);
X = vertcat(X, X_topNodes);

Geo.XgBottom = X_bottomIds;
Geo.XgTop = X_topIds;

xInternal = [1:nCells]';
Geo.nCells = length(xInternal);

%% Create tetrahedra
% [Twg_bottom] = CreateTetrahedra(trianglesConnectivity, neighboursNetwork, cellEdges, xInternal, X_bottomFaceIds, X_bottomVerticesIds);
% figure, tetramesh(Twg_bottom(any(ismember(Twg_bottom, xInternal), 2), :), X);
% 
% [Twg_top] = CreateTetrahedra(trianglesConnectivity_topoChanged, neighboursNetwork_topoChanged, cellEdges_topoChanged, xInternal, X_topFaceIds, X_topVerticesIds);
% Twg = vertcat(Twg_top, Twg_bottom);
% figure, tetramesh(Twg_top(any(ismember(Twg_top, xInternal), 2), :), X);
% figure, tetramesh(Twg(any(ismember(Twg, xInternal), 2), :), X);

X(X(:, 1) < 0 , 1) = 0;
X(X(:, 2) < 0 , 2) = 0;
X(X(:, 1) > 1 , 1) = 1;
X(X(:, 2) > 1 , 2) = 1;
tets = delaunayTriangulation(X);
Twg = tets.ConnectivityList;

%% Ghost cells and tets
Geo.XgID = setdiff(1:size(X, 1), xInternal);
Twg(all(ismember(Twg,Geo.XgID),2),:)=[];


%% After removing ghost tetrahedras, some nodes become disconnected,
% that is, not a part of any tetrahedra. Therefore, they should be
% removed from X
% Re-number the surviving tets
% uniqueTets = unique(Twg);
% Geo.XgID = Geo.nCells+1:length(uniqueTets);
% X    = X(uniqueTets,:);
% conv = zeros(size(X,1),1);
% conv(uniqueTets) = 1:size(X);
% Twg = conv(clTwg);

%% Build cells
[Geo] = BuildCells(Geo, Set, X, Twg);

%% Define upper and lower area threshold for remodelling
allFaces = [Geo.Cells.Faces];
allTris = [allFaces.Tris];
avgArea = mean([allTris.Area]);
stdArea = std([allTris.Area]);
Set.upperAreaThreshold = avgArea + stdArea;
Set.lowerAreaThreshold = avgArea - stdArea;

% TODO FIXME bad; PVM: better?
Geo.AssembleNodes = find(cellfun(@isempty, {Geo.Cells.AliveStatus})==0);
%% Define BarrierTri0
Set.BarrierTri0=realmax;
for c = 1:Geo.nCells
    Cell = Geo.Cells(c);
    for f = 1:length(Geo.Cells(c).Faces)
        Face = Cell.Faces(f);
        Set.BarrierTri0=min([vertcat(Face.Tris.Area); Set.BarrierTri0]);
    end
end
Set.BarrierTri0=Set.BarrierTri0/10;

end

