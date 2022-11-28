function [trianglesConnectivity, neighboursNetwork, cellEdges, verticesLocation, borderCells] = Build2DVoronoiFromImage(labelledImg, watershedImg, mainCells)
%BUILD3DVORONOITOPO Summary of this function goes here
%   Detailed explanation goes here

ratio = 5;
faceCentres = regionprops(labelledImg, 'centroid');
faceCentresVertices = fliplr(vertcat(faceCentres.Centroid));

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);

borderCellsAndMainCells = double(unique(vertcat(imgNeighbours{mainCells})));
borderCells = setdiff(borderCellsAndMainCells, mainCells);
uselessCells = setdiff(unique(labelledImg(:)), borderCellsAndMainCells);

labelledImg(ismember(labelledImg, uselessCells)) = 0;
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

[ verticesInfo ] = calculateVertices(labelledImg, watershedImg, imgNeighbours, 5);

%faceCentres = regionprops(labelledImg, 'centroid');
%faceCentresVertices = fliplr(vertcat(faceCentres.Centroid)) / imgSize;

totalCells = max(verticesInfo.connectedCells(:));
verticesInfo.PerCell = cell(totalCells, 1);

for numCell = borderCellsAndMainCells'
    verticesOfCell = find(any(ismember(verticesInfo.connectedCells, numCell), 2));
    verticesInfo.PerCell{numCell} = verticesOfCell;
    currentVertices = verticesInfo.location(verticesOfCell, :);
    currentConnectedCells = verticesInfo.connectedCells(verticesOfCell, :)';
    currentConnectedCells(currentConnectedCells == numCell) = [];
    currentConnectedCells = vertcat(currentConnectedCells(1:2:length(currentConnectedCells)), currentConnectedCells(2:2:length(currentConnectedCells)))';
    verticesInfo.edges{numCell, 1} = verticesOfCell(BoundaryOfCell(currentVertices, currentConnectedCells));
    if size(verticesInfo.edges{numCell, 1}, 1) < length(imgNeighbours{numCell})
        verticesInfo.edges{numCell, 1} = [];
    end
end

neighboursNetwork = [];

for numCell = borderCellsAndMainCells'
    currentNeighbours = double(imgNeighbours{numCell});
    currentCellNeighbours = [ones(length(currentNeighbours), 1) * numCell, currentNeighbours];
    
    neighboursNetwork = vertcat(neighboursNetwork, currentCellNeighbours);
end
trianglesConnectivity = double(verticesInfo.connectedCells);
cellEdges = verticesInfo.edges;
verticesLocation = verticesInfo.location;
end

