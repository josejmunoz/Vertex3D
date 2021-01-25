function [idCells] = findCentralCells(cellCentroids, numberOfCentralCells)
%FINDCENTRALCELLS Summary of this function goes here
%   Detailed explanation goes here

    centroidOfTissue = mean(cellCentroids, 1);
    [~, idCells] = pdist2(cellCentroids, centroidOfTissue, 'euclidean', 'Smallest', numberOfCentralCells);
end