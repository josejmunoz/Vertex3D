function Set = UpdateSet_F(Geo, Geo_0, Set)
%UPDATESET_F Summary of this function goes here
%   Detailed explanation goes here
    Set.f_Init = 0.75;
    Set.f = Set.f_Init * mean(cell2mat([Geo.Cells.Vol]) ./ cell2mat([Geo_0.Cells.Vol]))^3;
end

