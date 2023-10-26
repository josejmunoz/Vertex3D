function [quartets, percQuartets] = getFourFoldVertices(imgNeighbours)
%GETFOURFOLDVERTICES Summary of this function goes here
%   Detailed explanation goes here

quartets=buildQuartetsOfNeighs2D(imgNeighbours);
percQuartets=size(quartets,1)/length(imgNeighbours);

end