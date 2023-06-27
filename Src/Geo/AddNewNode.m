function [Geo, newNodeIDs] = AddNewNode(Geo, newXPositions, surroundingNodes)
%ADDNEWNODE Summary of this function goes here
%   Detailed explanation goes here
newNodeIDs = length(Geo.Cells)+1:length(Geo.Cells)+size(newXPositions, 1);
for newNodeID = 1:size(newXPositions, 1)
    Geo.Cells(newNodeIDs(newNodeID)).X = newXPositions(newNodeID, :);
    Geo.Cells(newNodeIDs(newNodeID)).T = [];
    Geo.Cells(newNodeIDs(newNodeID)).AliveStatus = [];
    Geo.Cells(newNodeIDs(newNodeID)).Area = [];
    Geo.Cells(newNodeIDs(newNodeID)).Area0 = [];
    Geo.Cells(newNodeIDs(newNodeID)).Vol = [];
    Geo.Cells(newNodeIDs(newNodeID)).Vol0 = [];
    Geo.Cells(newNodeIDs(newNodeID)).Y = [];
    Geo.Cells(newNodeIDs(newNodeID)).Faces = [];
    Geo.Cells(newNodeIDs(newNodeID)).cglobalIds = [];
    Geo.Cells(newNodeIDs(newNodeID)).globalIds = [];
    Geo.Cells(newNodeIDs(newNodeID)).ExternalLambda = [];
    Geo.Cells(newNodeIDs(newNodeID)).InternalLambda = [];
    Geo.Cells(newNodeIDs(newNodeID)).SubstrateLambda = [];
    Geo.XgID(end+1) = newNodeIDs(newNodeID);
    if DecideXgTopOrBottomByNeigh(Geo, surroundingNodes, newXPositions(newNodeID, :)) == 1
        Geo.XgTop(end+1) = newNodeIDs(newNodeID);
    else
        Geo.XgBottom(end+1) = newNodeIDs(newNodeID);
    end
end

end

