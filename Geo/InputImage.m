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

% Vertices of faces in depth: 
% x,y: mean([vertex1, vertex2], [], 2); z: 0
allEdges = vertcat(verticesInfo.edges{:});
allEdgesVertices = verticesInfo.location(allEdges, :);
allEdgesVertices = round(horzcat(mean(horzcat(allEdgesVertices(1:2:size(allEdgesVertices, 1), 1), allEdgesVertices(2:2:size(allEdgesVertices, 1), 1)), 2), mean(horzcat(allEdgesVertices(1:2:size(allEdgesVertices, 1), 2), allEdgesVertices(2:2:size(allEdgesVertices, 1), 2)), 2)));
Y = Y.Add([allEdgesVertices, zeros(size(allEdgesVertices, 1), 1)]);

%% Create Face Centres


%% Create nodes X from Y
% Cell nodes:
% x,y: vertices connected (edges); z: 0

% Ghost nodes:
% Above vertices (including faces) of top and bottom

%%%%% Go again to notes to see if anything is missing

end

