function [Geo, Set] = InitializeGeometry_3DVoronoi(Geo, Set)
%INITIALIZEGEOMETRY_3DVORONOI Summary of this function goes here
%   Detailed explanation goes here

nSeeds = Set.TotalCells + 30* Set.TotalCells;
imgDims = 3000;
lloydIterations = 100;
distorsion = 0;

rng default
x = randi(imgDims, nSeeds, 1);
y = randi(imgDims, nSeeds, 1);

seedsXY = horzcat(x,y);
seedsXY = unique(round(seedsXY, 2), 'rows');

%% Get central
distanceSeeds = pdist2(seedsXY, [imgDims/2 imgDims/2]);
[~, indices] = sort(distanceSeeds);
seedsXY = seedsXY(indices, :);

%% Homogeneize voronoi diagram
for numIter = 1:lloydIterations
    seedsXY(any(isnan(seedsXY), 2), :) = [];
    DT = delaunayTriangulation(seedsXY);
    [V, D] = voronoiDiagram(DT);
    for numCell = 1:size(seedsXY, 1)
        currentVertices = V(D{numCell}, :);
        seedsXY(numCell, :) = round(mean(currentVertices(all(~isinf(currentVertices), 2), :)));
    end
end

%% Get an image from it
img2D = zeros(imgDims, 'uint16');
for numCell = 1:size(seedsXY, 1)
    if all(seedsXY(numCell, :) > 0) && all(seedsXY(numCell, :) <= imgDims)
        img2D(seedsXY(numCell, 1), seedsXY(numCell, 2)) = 1;
    end
end

[distances, img2DLabelled] = bwdist(img2D);
watershedImg = watershed(distances, 8);

for numCell = 1:size(seedsXY, 1)
    if all(seedsXY(numCell, :) > 0) && all(seedsXY(numCell, :) <= imgDims)
        oldId = img2DLabelled(seedsXY(numCell, 1), seedsXY(numCell, 2));
        img2DLabelled(img2DLabelled == oldId) = numCell;
    end
end
img2DLabelled = uint16(img2DLabelled);

%% TODO: OBTAIN SHAPE OF THE CELLS TO ANALYSE THE ELLIPSE DIAMETER TO OBTAIN ITS REAL CELL HEIGHT
features2D = regionprops(img2DLabelled, 'all');
avgDiameter = mean([features2D(1:Set.TotalCells).MajorAxisLength]);
cellHeight = avgDiameter*Set.CellHeight;

%% TODO: Reorder here regarding the first cell (?)

%% Build 3D topology
[trianglesConnectivity, neighboursNetwork, cellEdges, verticesOfCell_pos, borderCells] = Build2DVoronoiFromImage(img2DLabelled, watershedImg, 1:Set.TotalCells);

seedsXY_topoChanged = [seedsXY(:, 1), seedsXY(:, 2) + rand(size(seedsXY, 1), 1)*distorsion];
% seedsXY_topoChanged(:, 2) = seedsXY_topoChanged(:, 2) - min(seedsXY_topoChanged(:, 2));
% seedsXY_topoChanged(:, 2) = seedsXY_topoChanged(:, 2) / max(seedsXY_topoChanged(:, 2));

[trianglesConnectivity_topoChanged, neighboursNetwork_topoChanged, cellEdges_topoChanged, verticesOfCell_pos_topoChanged] = Build2DVoronoiFromImage(img2DLabelled, watershedImg, 1:Set.TotalCells);

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

xInternal = [1:Set.TotalCells]';

%% Create tetrahedra
[Twg_bottom] = CreateTetrahedra(trianglesConnectivity, neighboursNetwork, cellEdges, xInternal, X_bottomFaceIds, X_bottomVerticesIds);
%[Twg_bottom, X, X_bottomIds] = upsampleTetMesh(Twg_bottom, X, X_bottomIds);
[Twg_top] = CreateTetrahedra(trianglesConnectivity_topoChanged, neighboursNetwork_topoChanged, cellEdges_topoChanged, xInternal, X_topFaceIds, X_topVerticesIds);
%[Twg_top, X, X_topIds] = upsampleTetMesh(Twg_top, X, X_topIds);
Twg = vertcat(Twg_top, Twg_bottom);

%% Fill Geo info
Geo.nCells = length(xInternal);

Geo.XgBottom = X_bottomIds;
Geo.XgTop = X_topIds;
Geo.XgLateral = setdiff(1:size(seedsXY, 1), xInternal);

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

%% Normalise Xs
X = X / imgDims;

%% Build cells
[Geo] = BuildCells(Geo, Set, X, Twg);

%% Define upper and lower area threshold for remodelling
allFaces = [Geo.Cells.Faces];
allTris = [allFaces.Tris];
avgArea = mean([allTris.Area]);
stdArea = std([allTris.Area]);
Set.upperAreaThreshold = avgArea + stdArea;
Set.lowerAreaThreshold = avgArea - stdArea;

%% Define border cells
Geo.BorderCells = borderCells;

Geo.BorderGhostNodes = setdiff(1:size(seedsXY, 1), 1:Geo.nCells);
Geo.BorderGhostNodes = [Geo.BorderGhostNodes'; setdiff(getNodeNeighbours(Geo, Geo.BorderGhostNodes), 1:Geo.nCells)];

% TODO FIXME bad; PVM: better?
Geo.AssembleNodes = find(cellfun(@isempty, {Geo.Cells.AliveStatus})==0);
%% Define BarrierTri0
Set.BarrierTri0=realmax;
Set.lmin0 = realmax;
edgeLengths_Top = [];
edgeLengths_Bottom = [];
edgeLengths_Lateral = [];
for c = 1:Geo.nCells
    Cell = Geo.Cells(c);
    for f = 1:length(Geo.Cells(c).Faces)
        Face = Cell.Faces(f);
        Set.BarrierTri0=min([vertcat(Face.Tris.Area); Set.BarrierTri0]);
        Set.lmin0=min([min(min(horzcat(vertcat(Face.Tris.LengthsToCentre), vertcat(Face.Tris.EdgeLength)))); Set.lmin0]);
        for tri = Face.Tris
            if tri.Location == 'Top'
                edgeLengths_Top(end+1) = ComputeEdgeLength(tri.Edge, Geo.Cells(c).Y);
            elseif tri.Location == 'Bottom'
                edgeLengths_Bottom(end+1) = ComputeEdgeLength(tri.Edge, Geo.Cells(c).Y);
            else
                edgeLengths_Lateral(end+1) = ComputeEdgeLength(tri.Edge, Geo.Cells(c).Y);
            end
        end
    end
    
    %Geo.Cells(c).Vol0 = mean([Geo.Cells(1:Geo.nCells).Vol]);
end
Geo.AvgEdgeLength_Top = mean(edgeLengths_Top);
Geo.AvgEdgeLength_Bottom = mean(edgeLengths_Bottom);
Geo.AvgEdgeLength_Lateral = mean(edgeLengths_Lateral);
Set.BarrierTri0=Set.BarrierTri0/5;
Set.lmin0 = Set.lmin0*10;

Geo.RemovedDebrisCells = [];

minZs = min(vertcat(Geo.Cells(1:Geo.nCells).Y));
Geo.CellHeightOriginal = abs(minZs(3));

end

