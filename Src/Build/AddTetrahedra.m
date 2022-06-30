function [Geo] = AddTetrahedra(Geo, newTets, Set)
%ADDTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here
for newTet = newTets'
    for numNode = newTet'
        if isempty(Geo.Cells(numNode).T) || all(~ismember(sort(Geo.Cells(numNode).T, 2), sort(newTet', 2), 'rows'))
            Geo.Cells(numNode).T(end+1, :) = newTet';
            if ~isempty(Geo.Cells(numNode).AliveStatus)
                Geo.Cells(numNode).Y(end+1, :) = ComputeY(vertcat(Geo.Cells(newTet).X), Geo.Cells(numNode).X, length([Geo.Cells(newTet).AliveStatus]) > 1, Set);
                Geo.numY = Geo.numY + 1;
            end
        end
    end
end
end

