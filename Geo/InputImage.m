function [] = InputImage(inputFile, cellHeight)
%INPUTIMAGE Summary of this function goes here
%   Detailed explanation goes here

img = imread(inputFile);

labelledImg = bwlabel(img, 8);

[imgNeighbours] = calculateNeighbours(labelledImg, 3);


end

