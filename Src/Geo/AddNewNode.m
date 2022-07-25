function [Geo, newNodeIDs] = AddNewNode(Geo, newXPosition)
%ADDNEWNODE Summary of this function goes here
%   Detailed explanation goes here
newNodeIDs = length(Geo.Cells)+1;
Geo.Cells(newNodeIDs).X = newXPosition;
Geo.XgID(end+1) = newNodeIDs;
%TODO: ADD ALSO TO BOTTOM OR TOP

end

