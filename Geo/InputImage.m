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
neighboursNetwork = [];

for numCell = 1:length(imgNeighbours)
    currentNeighbours = imgNeighbours{numCell};
    currentCellNeighbours = [ones(length(currentNeighbours), 1) * numCell, currentNeighbours];
    
    neighboursNetwork = vertcat(neighboursNetwork, currentCellNeighbours);
end
faceCentres = regionprops(labelledImg, 'centroid');
faceCentresVertices = fliplr(vertcat(faceCentres.Centroid)) / imgSize;

%% Create nodes X from Y
% Cell nodes:
X = horzcat(faceCentresVertices, zeros(size(faceCentresVertices, 1), 1));

% Ghost nodes:
% Above vertices (including faces) of top and bottom
XgTopFaceCentre = horzcat(faceCentresVertices, repmat(cellHeight, length(faceCentresVertices), 1));
XgBottomFaceCentre = horzcat(faceCentresVertices, repmat(-cellHeight, length(faceCentresVertices), 1));

X_bottomNodes = vertcat(XgBottomFaceCentre);
X_bottomIds = size(X, 1) + 1: size(X, 1) + size(X_bottomNodes, 1);
X = vertcat(X, X_bottomNodes);

X_topNodes = vertcat(XgTopFaceCentre);
X_topIds = size(X, 1) + 1: size(X, 1) + size(X_topNodes, 1);
X = vertcat(X, X_topNodes);

% Difference cell nodes and ghost nodes
xInternal = cellIdsAsInternal;
XgID = 1:size(X, 1);
XgID(xInternal) = [];

%% Create tetrahedra
% nodesDelaunay = delaunayTriangulation(X);
% Twg = nodesDelaunay.ConnectivityList;

nodesDelaunay = delaunayTriangulation(faceCentresVertices);
Twg_mid = nodesDelaunay.ConnectivityList;
Twg = [Twg_mid; Twg_mid + size(faceCentresVertices, 1); Twg_mid + (size(faceCentresVertices, 1)*2)];

%% Filtering of unnecessary nodes
% Remove Ghost tets 
Twg(sum(ismember(Twg, XgID), 2)~=2,:)=[];

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

