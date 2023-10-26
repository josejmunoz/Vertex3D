function [booleanIsCell] = IsCell(Geo, cell)
%ISCELL Summary of this function goes here
%   Detailed explanation goes here
    booleanIsCell = ~cellfun(@isempty, {Geo.Cells(cell).AliveStatus});
end

