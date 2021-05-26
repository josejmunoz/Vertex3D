function [quartets, percQuartets] = getFourFoldVertices(imgNeighbours, L_img)
%GETFOURFOLDVERTICES Summary of this function goes here
%   Detailed explanation goes here

totalCells=unique(L_img);
totalCells=totalCells(totalCells~=0);

try 
quartets=buildQuartetsOfNeighs2D(imgNeighbours);
% figure; imshow(ismember(L_img,quartets(:)))

catch
    quartets=[];
end
percQuartets=size(quartets,1)/length(imgNeighbours);

end

