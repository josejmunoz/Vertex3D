function [] = InputImage(inputFile, cellHeight)
%INPUTIMAGE Summary of this function goes here
%   Detailed explanation goes here

img = imread(inputFile);
labelledImg = bwlabel(1-img, 8);

ratio = 3;

[imgNeighbours] = calculateNeighbours(labelledImg, ratio);
[ verticesInfo ] = calculateVertices( labelledImg, imgNeighbours, ratio);

verticesInfo

end

