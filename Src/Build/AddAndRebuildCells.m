function [Geo_new] = AddAndRebuildCells(Geo, oldTets, newTets, Ynew, Set, updateMeasurements)
%ADDANDREBUILDCELLS Summary of this function goes here
%   Detailed explanation goes here
[Geo_new] = RemoveTetrahedra(Geo, oldTets);
[Geo_new] = AddTetrahedra(Geo_new, Geo, newTets, Ynew, Set);
Geo_new = Rebuild(Geo_new, Set);
Geo_new   = BuildGlobalIds(Geo_new);
%if updateMeasurements
    Geo_new   = UpdateMeasures(Geo_new);
%end

Geo_new = CheckYsAndFacesHaveNotChanged(Geo, newTets, Geo_new);
end

