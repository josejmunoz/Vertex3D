function [Ynew] = RecalculateYs(Geo, Tnew, mainNodesToConnect, Set)
%RECALCULATEYS Summary of this function goes here
%   Detailed explanation goes here

allTs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).T);
allYs = vertcat(Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).Y);
Ynew = [];
for numTet = 1:size(Tnew, 1)
    mainNode_current = mainNodesToConnect(ismember(mainNodesToConnect, Tnew(numTet, :)));
    Ynew(end+1, :) = ComputeY(Geo, Tnew(numTet, :), Geo.Cells(mainNode_current(1)).X, Set);
end
end

