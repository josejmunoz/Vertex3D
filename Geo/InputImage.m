function [X, Y0, Y,T, Tetrahedra_weights,XgID,Cell,Cn,Cv,Yn,SCn,Set] = InputImage(Set)
%INPUTIMAGE Summary of this function goes here
%   Detailed explanation goes here

%% Processed input image to obtain the information
img = imread(Set.InputSegmentedImage);
img(1:end, 1) = 1;
img(1:end, end) = 1;
img(1, 1:end) = 1;
img(end, 1:end) = 1;
labelledImg = bwlabel(1-img, 8);

imgSize = size(labelledImg, 1);
totalCells = Set.TotalCells;

ratio = 5;
faceCentres = regionprops(labelledImg, 'centroid');
faceCentresVertices = fliplr(vertcat(faceCentres.Centroid));
cellIdsAsInternal = findCentralCells(faceCentresVertices, size(faceCentresVertices, 1));

newLabelledImg = zeros(size(labelledImg));
for numCell = 1:size(faceCentresVertices, 1)
    %newLabelledImg(ismember(labelledImg, numCell)) = cellIdsAsInternal(numCell);
    newLabelledImg(ismember(labelledImg, cellIdsAsInternal(numCell))) = numCell;
end

labelledImg = newLabelledImg;

cellIdsAsInternal = findCentralCells(faceCentresVertices, totalCells);


cellArea = regionprops(labelledImg, 'Area');
avgCellArea = mean(vertcat(cellArea(1:totalCells).Area))/(imgSize^2);
cellHeight = Set.CellHeight * avgCellArea;

cellIdsAsInternal = 1:totalCells;

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);

%% Remove quartets
[quartets] = getFourFoldVertices(imgNeighbours, labelledImg);
for numQuartets = 1:size(quartets, 1)
    currentCentroids = faceCentresVertices(quartets(numQuartets, :), :);
    distanceBetweenCentroids = squareform(pdist(currentCentroids));
    [maxDistance] = max(distanceBetweenCentroids(:));
    [row, col] = find(distanceBetweenCentroids == maxDistance);
    
    % Remove first neighbour from the furthest pair of neighbour
    currentNeighs = imgNeighbours{quartets(numQuartets, col(1))};
    currentNeighs(currentNeighs == quartets(numQuartets, row(1))) = [];
    imgNeighbours{quartets(numQuartets, col(1))} = currentNeighs;
    
    % Remove the second of the same pair
    currentNeighs = imgNeighbours{quartets(numQuartets, row(1))};
    currentNeighs(currentNeighs == quartets(numQuartets, col(1))) = [];
    imgNeighbours{quartets(numQuartets, row(1))} = currentNeighs;
end

[ verticesInfo ] = calculateVertices( labelledImg, imgNeighbours, ratio);
%% Slightly change the positions of repated vertices (or fourfold verties)
[~, ia] = unique(verticesInfo.location, 'rows');
repeatedVertices = find(ismember(1:size(verticesInfo.location, 1), ia) == 0);
verticesInfo.location(repeatedVertices, 1) = verticesInfo.location(repeatedVertices, 1)+1;

vertex2D(ismember(1:size(vertex2D, 1), ia) == 0, :);


faceCentres = regionprops(labelledImg, 'centroid');
faceCentresVertices = fliplr(vertcat(faceCentres.Centroid)) / imgSize;

totalCells = max(verticesInfo.connectedCells(:));
verticesInfo.PerCell = cell(totalCells, 1);

for numCell = 1:totalCells
    verticesOfCell = find(any(ismember(verticesInfo.connectedCells, numCell), 2));
    verticesInfo.PerCell{numCell} = verticesOfCell;
    currentVertices = verticesInfo.location(verticesOfCell, :);
    currentConnectedCells = verticesInfo.connectedCells(verticesOfCell, :)';
    currentConnectedCells(currentConnectedCells == numCell) = [];
    currentConnectedCells = vertcat(currentConnectedCells(1:2:length(currentConnectedCells)), currentConnectedCells(2:2:length(currentConnectedCells)))';
    verticesInfo.edges{numCell, 1} = verticesOfCell(boundaryOfCell(currentVertices, currentConnectedCells));
    if size(verticesInfo.edges{numCell, 1}, 1) < length(imgNeighbours{numCell})
        verticesInfo.edges{numCell, 1} = [];
    end
end

neighboursNetwork = [];

for numCell = 1:length(imgNeighbours)
    currentNeighbours = imgNeighbours{numCell};
    currentCellNeighbours = [ones(length(currentNeighbours), 1) * numCell, currentNeighbours];
    
    neighboursNetwork = vertcat(neighboursNetwork, currentCellNeighbours);
end

% figure, imshow(labelledImg, colorcube(600))
% hold on;
% for numVertex = 1:size(verticesInfo.location, 1)
%     %plot(round(faceCentresVertices(numVertex, 2)), round(faceCentresVertices(numVertex, 1)), 'bo');
%     plot(round(verticesInfo.location(numVertex, 2)), round(verticesInfo.location(numVertex, 1)), 'rx');
%     hold on;
% end
% figure, imshow(ismember(labelledImg, find(cellfun(@isempty, verticesInfo.edges) == 0)).*labelledImg, colorcube(600))

%% Create nodes X from Y
% Cell nodes:
% x,y: centroidsOfCells; z: 0
nonEmptyCells = cellfun(@isempty, verticesInfo.edges) == 0;

nonEmptyCells = zeros(size(nonEmptyCells)) == 1;
nonEmptyCells(cellIdsAsInternal) = 1;

X = horzcat(faceCentresVertices, zeros(size(nonEmptyCells)));

vertex2D = verticesInfo.location / imgSize;

% Ghost nodes:
% Above vertices (including faces) of top and bottom
XgTopFaceCentre = horzcat(faceCentresVertices, repmat(cellHeight, length(faceCentresVertices), 1));
XgBottomFaceCentre = horzcat(faceCentresVertices, repmat(-cellHeight, length(faceCentresVertices), 1));
XgTopVertices = [vertex2D, repmat(cellHeight, size(vertex2D, 1), 1)];
XgBottomVertices = [vertex2D, repmat(-cellHeight, size(vertex2D, 1), 1)];

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


% Difference cell nodes and ghost nodes
xInternal = find(nonEmptyCells);
XgID = horzcat(find(nonEmptyCells == 0)', (length(faceCentresVertices)+1):size(X, 1));

%% Create tetrahedra
trianglesConnectivity = verticesInfo.connectedCells;
[Twg_bottom] = createTetrahedra(trianglesConnectivity, neighboursNetwork, verticesInfo.edges, xInternal, X_bottomFaceIds, X_bottomVerticesIds);
[Twg_top] = createTetrahedra(trianglesConnectivity, neighboursNetwork, verticesInfo.edges, xInternal, X_topFaceIds, X_topVerticesIds);

Twg = vertcat(Twg_top, Twg_bottom);

%% Filtering of unnecessary nodes
% Remove Ghost tets 
Twg(all(ismember(Twg,XgID),2),:)=[];

%Centre Nodal position at (0,0)
X(:,1)=X(:,1)-mean(X(:,1));
X(:,2)=X(:,2)-mean(X(:,2));
X(:,3)=X(:,3)-mean(X(:,3));

%% Check order of tetrahedrons
[Twg] = CheckTetrahedronOrder(Twg, X);

%% Create vertices Y
%Y_new=GetYFromX(X,XgID,Twg,cellHeight/2);
Y_new = zeros(size(Twg, 1), 3);
Tetrahedra_weights = ones(size(Twg, 1), 3);
for numTetrahedron = 1:size(Twg, 1)
    Y_new(numTetrahedron, :) = mean(X(Twg(numTetrahedron, :), :));
    if Y_new(numTetrahedron, 3) > 0
        Tetrahedra_weights(numTetrahedron, :) = [1 1 Y_new(numTetrahedron, 3)/(cellHeight/2)];
    elseif Y_new(numTetrahedron, 3) < 0
        Tetrahedra_weights(numTetrahedron, :) = [1 1 Y_new(numTetrahedron, 3)/(-cellHeight/2)];
    end
    Y_new(numTetrahedron, :) = Y_new(numTetrahedron, :) ./ Tetrahedra_weights(numTetrahedron, :);
end

Y=DynamicArray(ceil(size(Y_new,1)*1.5),size(Y_new,2));
Y=Y.Add(Y_new);

% % sum(any(ismember(Twg, xInternal(2)), 2))
% index = any(ismember(Twg, xInternal), 2);
% %figure, tetramesh(Twg(index, :), X);
% figure, plot3(Y.DataRow(index, 1), Y.DataRow(index, 2), Y.DataRow(index, 3), 'rx');
% %indexVertices = sum(ismember(Twg, borderPairs), 2) == 1;
% indexVertices = sum(ismember(Twg, borderPairs(1:22, 2)), 2) >= 1;
% hold on, plot3(Y.DataRow(indexVertices, 1), Y.DataRow(indexVertices, 2), Y.DataRow(indexVertices, 3), 'bo');


%% Create cells
xInternal = xInternal';
[Cv,Cell]=BuildCells(Twg,Y,X,xInternal, cellHeight, false);

borderPairs = neighboursNetwork(sum(ismember(neighboursNetwork, xInternal), 2) == 1, :);
borderPairs = unique(sort(borderPairs, 2), 'rows');
Cell.BorderVertices = find(sum(ismember(Twg, borderPairs(:, 2)), 2) >= 1);
% Add facecentres
Cell.BorderVertices = [Cell.BorderVertices; -find(ismember(Cell.AllFaces.Nodes, borderPairs, 'rows'))];
Cell.BorderCells = ismember(Cell.Int, borderPairs(:));
Set.NumMainV=Y.n;
Set.NumAuxV=Cell.FaceCentres.n;
Set.NumCellCentroid = Cell.n;
Set.NumTotalV=Set.NumMainV + Set.NumAuxV + Set.NumCellCentroid;
Set.NumXs = size(X, 1);
Cn=BuildCn(Twg);

Cell.AllFaces=Cell.AllFaces.ComputeAreaTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.ComputePerimeterTri(Y.DataRow,Cell.FaceCentres.DataRow);
Cell.AllFaces=Cell.AllFaces.CheckInteriorFaces(XgID);

Yn = Y;
Y0 = Y;
SCn = Cell.FaceCentres;
Cell.Centre_n = Cell.Centre;

T=DynamicArray(ceil(size(Twg,1)*1.5),size(Twg,2));
T=T.Add(Twg);

%[X]=GetXFromY(Cell,X,T,Y,XgID,Set);

% figure,
% tetramesh(T.DataRow(any(ismember(T.DataRow, 1), 2), :), XNew);


% Set.BarrierTri0=realmax; 
% for i=1:Cell.n
%     Set.BarrierTri0=min([Cell.SAreaTri{i}; Set.BarrierTri0]);
% end
% Set.BarrierTri0=Set.BarrierTri0/10;

allPerimeters = vertcat(Cell.AllFaces.PerimeterTri{:});
allAreas = vertcat(Cell.AllFaces.AreaTri{:});
Set.BarrierTri0=min(allPerimeters .* allAreas)/2;

[Cell,Y]=CheckOrderingOfTriangulaiton(Cell,Y,Set);

end

