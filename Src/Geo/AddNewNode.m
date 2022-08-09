function [Geo, newNodeIDs] = AddNewNode(Geo, newXPositions)
%ADDNEWNODE Summary of this function goes here
%   Detailed explanation goes here
newNodeIDs = length(Geo.Cells)+1:length(Geo.Cells)+size(newXPositions, 1);
for newNodeID = 1:size(newXPositions, 1)
    Geo.Cells(newNodeIDs(newNodeID)).X = newXPositions(newNodeID, :);
    Geo.XgID(end+1) = newNodeIDs(newNodeID);
    %TODO: ADD ALSO TO BOTTOM OR TOP
end

end

