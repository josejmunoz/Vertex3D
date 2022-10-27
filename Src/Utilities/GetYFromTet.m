function [Y, YId] = GetYFromTet(Geo, tet)
%GETYFROMTET Summary of this function goes here
%   Detailed explanation goes here
    mainNodes = tet(~cellfun(@isempty, {Geo.Cells(tet).AliveStatus}));
    Ts = Geo.Cells(mainNodes(1)).T;
    [foundYs]= ismember(sort(Ts, 2), sort(tet, 2), 'rows');
    Y = Geo.Cells(mainNodes(1)).Y(foundYs, :);
    YId = find(foundYs);
end

