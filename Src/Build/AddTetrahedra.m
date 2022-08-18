function [Geo] = AddTetrahedra(Geo, newTets, Set)
%ADDTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here
for newTet = newTets'
    for numNode = newTet'
        if isempty(Geo.Cells(numNode).T) || all(~ismember(sort(Geo.Cells(numNode).T, 2), sort(newTet', 2), 'rows'))
            Geo.Cells(numNode).T(end+1, :) = newTet';
            if ~isempty(Geo.Cells(numNode).AliveStatus) && exist('Set', 'var')
                if isequal(Set.InputGeo, 'Voronoi')
                    Geo.Cells(numNode).Y(end+1, :) = ComputeY(vertcat(Geo.Cells(newTet).X), Geo.Cells(numNode).X, 0, Set);
                    if sum(ismember(Geo.Cells(numNode).T(end, :), Geo.XgTop)) > 0
                        Geo.Cells(numNode).Y(end, 3) = Geo.Cells(numNode).Y(end, 3) / (sum(ismember(newTet, Geo.XgTop))/2);
                    elseif sum(ismember(Geo.Cells(numNode).T(end, :), Geo.XgBottom)) > 0
                        Geo.Cells(numNode).Y(end, 3) = Geo.Cells(numNode).Y(end, 3) / (sum(ismember(newTet, Geo.XgBottom))/2);
                    end
                else
                    Geo.Cells(numNode).Y(end+1, :) = ComputeY(vertcat(Geo.Cells(newTet).X), Geo.Cells(numNode).X, length([Geo.Cells(newTet).AliveStatus]), Set);
                end
                
                if any(isinf(Geo.Cells(numNode).Y(end, :)))
                    disp('mierda');
                end
                Geo.numY = Geo.numY + 1;
            end
        end
    end   
end
end

