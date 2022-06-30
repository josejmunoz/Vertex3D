function [convexTet] = CheckConvexityTets(nodeTets, newTets, Geo, nodeItShouldBeConnected)
%CHECKCONVEXITYTETS Summary of this function goes here
%   Detailed explanation goes here

%Combine the two new tets created
uniqueNewAssigned = unique([nodeTets; newTets]);
%Obtain the convex shape connecting the tets
shp_newAssigned = uniqueNewAssigned(boundary(vertcat(Geo.Cells(uniqueNewAssigned).X), 0));
%Obtain how they are connected
neighbours_Assigned = unique(shp_newAssigned(any(ismember(shp_newAssigned, nodeItShouldBeConnected), 2), :));

%It would be convex if the connections correspond to their assigned node
%(nodeTets)
convexTet = all(ismember(nodeTets, neighbours_Assigned));
end

