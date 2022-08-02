function [] = InitializeGeometry_3DVoronoi()
%INITIALIZEGEOMETRY_3DVORONOI Summary of this function goes here
%   Detailed explanation goes here

nSeeds = 100;
distorsion = 1;
cellHeight = 1;

rng default
x = rand(nSeeds, 1);
y = rand(nSeeds, 1);

seedsXY = horzcat(x,y);
%% Get central
distanceSeeds = pdist2(seedsXY, [0.5 0.5]);
[~, indices] = sort(distanceSeeds);
seedsXY = seedsXY(indices, :);

[trianglesConnectivity, neighboursNetwork, cellEdges, verticesOfCell_pos] = Build3DVoronoiTopo(seedsXY);

seedsXY_topoChanged = [seedsXY(:, 1), seedsXY(:, 2) + distorsion];
seedsXY_topoChanged(:, 2) = seedsXY_topoChanged(:, 2) - min(seedsXY_topoChanged(:, 2));
seedsXY_topoChanged(:, 2) = seedsXY_topoChanged(:, 2) / max(seedsXY_topoChanged(:, 2));

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

X_IDs.bottomVerticesIds = X_bottomVerticesIds;
X_IDs.bottomFaceIds = X_bottomFaceIds;
X_IDs.topVerticesIds = X_topVerticesIds;
X_IDs.topFaceIds = X_topFaceIds;

xInternal = [1:20]';

%% Create tetrahedra
X(X(:, 1) < 0 , 1) = 0;
X(X(:, 2) < 0 , 2) = 0;
X(X(:, 1) > 1 , 1) = 1;
X(X(:, 2) > 1 , 2) = 1;
tets = delaunayTriangulation(X);
Twg = tets.ConnectivityList;

figure, tetramesh(Twg, X);
figure, tetramesh(Twg(any(ismember(Twg, xInternal), 2), :), X);
Twg

end

