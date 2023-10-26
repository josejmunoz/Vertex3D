function [sphericity] = ComputeCellSphericity(Cell)
%COMPUTECELLSPHERICITY Summary of this function goes here
%   Detailed explanation goes here
% Based on https://sciencing.com/calculate-sphericity-5143572.html
sphericity = (pi^(1/3) * (6 * Cell.Vol) ^ (2/3)) / Cell.Area;
end

