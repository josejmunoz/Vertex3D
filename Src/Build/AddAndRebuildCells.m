function [Geo_new] = AddAndRebuildCells(Geo, oldTets, newTets, Ynew, Set, updateMeasurements)
%ADDANDREBUILDCELLS Summary of this function goes here
%   Detailed explanation goes here
[Geo_new] = RemoveTetrahedra(Geo, oldTets);
[Geo_new] = AddTetrahedra(Geo_new, Geo, newTets, Ynew, Set);
Geo_new = Rebuild(Geo_new, Set);
Geo_new   = BuildGlobalIds(Geo_new);

Geo_new = CheckYsAndFacesHaveNotChanged(Geo, newTets, Geo_new);

%if updateMeasurements
    Geo_new   = UpdateMeasures(Geo_new);
%end

%% Check here how many neighbours they're loosing and winning and change the number of lambdaA_perc accordingly
neighbours_init = [];
for cell = Geo.Cells(1:Geo.nCells)
    neighbours_init(end+1) = length(getNodeNeighbours(Geo, cell.ID));
end

neighbours_end = [];
for cell = Geo_new.Cells(1:Geo_new.nCells)
    neighbours_end(end+1) = length(getNodeNeighbours(Geo_new, cell.ID));
end

difference = neighbours_init - neighbours_end;

for numCell = 1:Geo.nCells
    Geo.Cells(numCell).lambdaB_perc = Geo.Cells(numCell).lambdaB_perc - (0.01 * difference(numCell));
end

