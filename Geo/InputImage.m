function [] = InputImage(inputFile, cellHeight)
%INPUTIMAGE Summary of this function goes here
%   Detailed explanation goes here

%% Processed input image to obtain the information
img = imread(inputFile);
img(1:end, 1) = 1;
img(1:end, end) = 1;
img(1, 1:end) = 1;
img(end, 1:end) = 1;
labelledImg = bwlabel(1-img, 8);

ratio = 5;

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);
[ verticesInfo ] = calculateVertices( labelledImg, imgNeighbours, ratio);

faceCentres = regionprops(labelledImg, 'centroid');

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

% figure, imshow(labelledImg, colorcube(600))
% figure, imshow(ismember(labelledImg, find(cellfun(@isempty, verticesInfo.edges) == 0)).*labelledImg, colorcube(600))

%% Create vertices Y
Y=DynamicArray(size(verticesInfo.location, 1)*3, 3);

% Top vertices:
% x,y: imgNeighbours; z:cellHeight/2
vertex2D = verticesInfo.location;
Y = Y.Add([vertex2D, repmat(cellHeight/2, size(vertex2D, 1), 1)]);

% Bottom vertices:
% x,y: imgNeighbours; z:-cellHeight/2
Y = Y.Add([vertex2D, repmat(-cellHeight/2, size(vertex2D, 1), 1)]);

% % Vertices of faces in depth:
% % x,y: mean([vertex1, vertex2], [], 2); z: 0
% allEdges = vertcat(verticesInfo.edges{:});
% allEdgesVertices = verticesInfo.location(allEdges, :);
% allEdgesVertices = round(horzcat(mean(horzcat(allEdgesVertices(1:2:size(allEdgesVertices, 1), 1), allEdgesVertices(2:2:size(allEdgesVertices, 1), 1)), 2), mean(horzcat(allEdgesVertices(1:2:size(allEdgesVertices, 1), 2), allEdgesVertices(2:2:size(allEdgesVertices, 1), 2)), 2)));
% Y = Y.Add([allEdgesVertices, zeros(size(allEdgesVertices, 1), 1)]);


%% Create nodes X from Y
% Cell nodes:
% x,y: centroidsOfCells; z: 0
nonEmptyCells = cellfun(@isempty, verticesInfo.edges) == 0;
totalRegularCells = sum(nonEmptyCells);
xInternal = 1:totalRegularCells;
X = horzcat(vertcat(faceCentres(nonEmptyCells).Centroid), zeros(totalRegularCells, 1));

% Ghost nodes:
% Above vertices (including faces) of top and bottom
XgTopFaceCentre = horzcat(vertcat(faceCentres(nonEmptyCells).Centroid), repmat(cellHeight, size(vertex2D, 1), 1));
XgBottomFaceCentre = horzcat(vertcat(faceCentres(nonEmptyCells).Centroid), repmat(-cellHeight, size(vertex2D, 1), 1));
XgTopVertices = [vertex2D, repmat(cellHeight, size(vertex2D, 1), 1)];
XgBottomVertices = [vertex2D, repmat(-cellHeight, size(vertex2D, 1), 1)];


XgID = totalRegularCells+1:size(X, 1);

%% Create tetrahedra
Twg=delaunay(X);

% Remove Ghost tets 
Twg(all(ismember(Twg,XgID),2),:)=[];

% Remove addition nodes- not used
newTwg=zeros(size(Twg));
newX=zeros(size(X));
newXgID=zeros(size(X,1),1);
aux1=zeros(size(X,1),1);
aux2=1;
aux3=1;
for i=1:size(X,1)
    if ismember(i,Twg)
       newTwg(ismember(Twg,i))=aux2;
       newX(aux2,:)=X(i,:);
       aux1(i)=aux2;
       if ismember(i,XgID)
           newXgID(aux3)=aux2;
           aux3=aux3+1;
       end 
       aux2=aux2+1;
    end 
end 
X=newX(1:aux2-1,:);
XgID=newXgID(1:aux3-1);
Twg=newTwg;
% Keep only nodes-ghost, the others are already from verticesInfo

%% Create cells
[Cv,Cell,SharedFaces]=BuildCells(Twg,Y,X,xInternal,H);

%%%%% Go again to notes to see if anything is missing

end

