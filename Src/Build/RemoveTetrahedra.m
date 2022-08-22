function [Geo, oldYs] = RemoveTetrahedra(Geo, removingTets)
%REMOVETETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here
oldYs = [];
for removingTet = removingTets'
    for numNode = removingTet'
        idToRemove = ismember(sort(Geo.Cells(numNode).T, 2), sort(removingTet', 2), 'rows');
        Geo.Cells(numNode).T(idToRemove, :) = [];
        if ~isempty(Geo.Cells(numNode).AliveStatus)
            oldYs = [oldYs, Geo.Cells(numNode).Y(idToRemove, :)];
            Geo.Cells(numNode).Y(idToRemove, :) = [];
            Geo.numY = Geo.numY - 1;
        end
    end
end
end

