function [Geo] = AddTetrahedra(Geo, newTets, Set)
%ADDTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here
for newTet = newTets'
    for numNode = newTet'
        if isempty(Geo.Cells(numNode).T) || all(~ismember(sort(Geo.Cells(numNode).T, 2), sort(newTet', 2), 'rows'))
            Geo.Cells(numNode).T(end+1, :) = newTet';
            if ~isempty(Geo.Cells(numNode).AliveStatus) && exist('Set', 'var')
                Geo.Cells(numNode).Y(end+1, :) = ComputeY(vertcat(Geo.Cells(newTet).X), Geo.Cells(numNode).X, length([Geo.Cells(newTet).AliveStatus]), Set);
                Geo.numY = Geo.numY + 1;
                
                if isequal(Set.InputGeo, 'Voronoi')
                    if sum(ismember(Geo.Cells(numNode).T, Geo.XgTop)) > 0
                        Geo.Cells(numNode).Y(end, 3) = Geo.Cells(numNode).Y(end, 3) / (sum(ismember(newTet, Geo.XgTop))/2);
                    elseif sum(ismember(Geo.Cells(numNode).T, Geo.XgBottom)) > 0
                        Geo.Cells(numNode).Y(end, 3) = Geo.Cells(numNode).Y(end, 3) / (sum(ismember(newTet, Geo.XgBottom))/2);
                    end
                end
            end
        end
    end
end
end

