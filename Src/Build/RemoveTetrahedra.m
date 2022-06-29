function [Geo] = RemoveTetrahedra(Geo, removingTets)
%REMOVETETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here
for removingTet = removingTets'
    for numNode = removingTet'
        idToRemove = ismember(sort(Geo.Cells(numNode).T, 2), sort(removingTet', 2), 'rows');
        Geo.Cells(numNode).T(idToRemove, :) = [];
        if ~isempty(Geo.Cells(numNode).AliveStatus)
            Geo.Cells(numNode).Y(idToRemove, :) = [];
        end
    end
end
end

