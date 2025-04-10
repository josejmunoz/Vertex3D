function [trianglesConnectivity, neighboursNetwork, cellEdges, verticesLocation, borderCells, borderOfborderCellsAndMainCells] = Build2DVoronoiFromImage(labelledImg, watershedImg, mainCells)
%BUILD3DVORONOITOPO Summary of this function goes here
%   Detailed explanation goes here

ratio = 2;

labelledImg(watershedImg == 0) = 0;

% Create a mask for the edges with ID 0
edgeMask = labelledImg == 0;

% Get the closest labeled polygon for each edge pixel
closestID = imdilate(labelledImg, true(ratio));

filledImage = closestID;
filledImage(~edgeMask) = labelledImg(~edgeMask);

labelledImg = filledImage;

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);

borderCellsAndMainCells = double(unique(vertcat(imgNeighbours{mainCells})));
borderGhostCells = setdiff(borderCellsAndMainCells, mainCells); 
borderCells = intersect(mainCells, double(unique(vertcat(imgNeighbours{borderGhostCells}))));

borderOfborderCellsAndMainCells = double(unique(vertcat(imgNeighbours{borderCellsAndMainCells})))';
labelledImg(~ismember(labelledImg, 1:max(borderOfborderCellsAndMainCells))) = 0;
[imgNeighbours] = calculateNeighbours(labelledImg, ratio);

% %% Remove quartets
[quartets] = getFourFoldVertices(imgNeighbours);
faceCentres = regionprops(labelledImg, 'centroid');
faceCentresVertices = fliplr(vertcat(faceCentres.Centroid));
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

[ verticesInfo ] = calculateVertices(labelledImg, imgNeighbours, ratio);

%faceCentres = regionprops(labelledImg, 'centroid');
%faceCentresVertices = fliplr(vertcat(faceCentres.Centroid)) / imgSize;

totalCells = max(borderCellsAndMainCells);
verticesInfo.PerCell = cell(totalCells, 1);

for numCell = 1:max(mainCells)
    verticesOfCell = find(any(ismember(verticesInfo.connectedCells, numCell), 2));
    verticesInfo.PerCell{numCell} = verticesOfCell;
    currentVertices = verticesInfo.location(verticesOfCell, :);
    currentConnectedCells = verticesInfo.connectedCells(verticesOfCell, :)';
    currentConnectedCells(currentConnectedCells == numCell) = [];
    currentConnectedCells = vertcat(currentConnectedCells(1:2:length(currentConnectedCells)), currentConnectedCells(2:2:length(currentConnectedCells)))';
    verticesInfo.edges{numCell, 1} = verticesOfCell(BoundaryOfCell(currentVertices, currentConnectedCells));
    assert(size(verticesInfo.edges{numCell, 1}, 1) == length(imgNeighbours{numCell}), strcat('Error missing vertices of neighbours in ', num2str(numCell)))
end

neighboursNetwork = [];

for numCell = 1:max(mainCells)
    currentNeighbours = double(imgNeighbours{numCell});
    currentCellNeighbours = [ones(length(currentNeighbours), 1) * numCell, currentNeighbours];
    
    neighboursNetwork = vertcat(neighboursNetwork, currentCellNeighbours);
end

%% Final assigning
trianglesConnectivity = double(verticesInfo.connectedCells);
cellEdges = verticesInfo.edges;
verticesLocation = verticesInfo.location;
end

