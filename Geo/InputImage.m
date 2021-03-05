function [] = InputImage(inputFile, cellHeight)
%INPUTIMAGE Summary of this function goes here
%   Detailed explanation goes here

%% Processed input image to obtain the information
img = imread(inputFile);
labelledImg = bwlabel(1-img, 8);

ratio = 3;

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);
[ verticesInfo ] = calculateVertices( labelledImg, imgNeighbours, ratio);

faceCentres = regionprops(labelledImg, 'centroid');

totalCells = max(verticesInfo.connectedCells(:));
verticesInfo.PerCell = cell(totalCells, 1);

for numCell = 2%1:totalCells
    verticesOfCell = find(any(ismember(verticesInfo.connectedCells, numCell), 2));
    
end

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


%% Create Face Centres


%% Create nodes X from Y
% Cell nodes:
% x,y: vertices connected (edges); z: 0

% Ghost nodes:
% Above vertices (including faces) of top and bottom

%%%%% Go again to notes to see if anything is missing

end

