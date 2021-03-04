function [imgNeighbours] = calculateNeighbours(labelledImg, ratioStrel)
%CALCULATENEIGHBOURS Summary of this function goes here
%   Detailed explanation goes here

    se = strel('disk',ratioStrel);
    
    cells=sort(unique(labelledImg));
    cells=cells(2:end);                  %% Deleting cell 0 from range
    imgNeighbours=cell(length(cells),1);

    for cel = 1:length(cells)
        BW = bwperim(labelledImg==cells(cel));
        BW_dilate=imdilate(BW,se);
        neighs=unique(labelledImg(BW_dilate==1));
        imgNeighbours{cells(cel)}=neighs((neighs ~= 0 & neighs ~= cells(cel)));
    end
end

