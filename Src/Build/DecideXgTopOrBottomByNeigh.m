function [location] = DecideXgTopOrBottomByNeigh(Geo, surroundingNodes, newNodePosition)
%DECIDEXGTOPORBOTTOMBYNEIGH Summary of this function goes here
%   Detailed explanation goes here

percentageOfBottom = sum(ismember(surroundingNodes, Geo.XgBottom))/numel(surroundingNodes);
percentageOfTop = sum(ismember(surroundingNodes, Geo.XgTop))/numel(surroundingNodes);

if percentageOfBottom > 2*percentageOfTop
    location = 2; %% Bottom
elseif 2*percentageOfBottom < percentageOfTop
    location = 1; %% Top
else %% Need further tests
    [~, idClosest] = pdist2(newNodePosition, vertcat(Geo.Cells(surroundingNodes).X), 'euclidean', 'smallest', 1);
    if ismember(surroundingNodes(idClosest), Geo.XgTop)
        location = 1; %% Top
    else
        location = 2; %% Bottom
    end
end

