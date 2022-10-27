function [Geo, newNodeIDs] = AddNewNode(Geo, newXPositions, surroundingNodes)
%ADDNEWNODE Summary of this function goes here
%   Detailed explanation goes here
newNodeIDs = length(Geo.Cells)+1:length(Geo.Cells)+size(newXPositions, 1);
for newNodeID = 1:size(newXPositions, 1)
    Geo.Cells(newNodeIDs(newNodeID)).X = newXPositions(newNodeID, :);
    Geo.XgID(end+1) = newNodeIDs(newNodeID);
    if DecideXgTopOrBottomByNeigh(Geo, surroundingNodes, newXPositions(newNodeID, :)) == 1
        Geo.XgTop(end+1) = newNodeIDs(newNodeID);
    else
        Geo.XgBottom(end+1) = newNodeIDs(newNodeID);
    end
end

end

