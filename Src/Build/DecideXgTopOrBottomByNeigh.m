function [location] = DecideXgTopOrBottomByNeigh(Geo, surroundingNodes, newNodePosition)
%DECIDEXGTOPORBOTTOMBYNEIGH Summary of this function goes here
%   Detailed explanation goes here

bottomGhostNodes = ismember(surroundingNodes, Geo.XgBottom);
topGhostNodes = ismember(surroundingNodes, Geo.XgTop);
percentageOfBottom = sum(bottomGhostNodes)/numel(surroundingNodes);
percentageOfTop = sum(topGhostNodes)/numel(surroundingNodes);

if percentageOfBottom > 2*percentageOfTop
    location = 2; %% Bottom
elseif 2*percentageOfBottom < percentageOfTop
    location = 1; %% Top
else %% Need further tests
    [distances] = pdist2(vertcat(Geo.Cells(surroundingNodes).X), newNodePosition, 'euclidean');
    if mean(distances(topGhostNodes)) < mean(distances(bottomGhostNodes))
        location = 1; %% Top
    else
        location = 2; %% Bottom
    end
end

