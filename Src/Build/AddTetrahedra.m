function [Geo] = AddTetrahedra(Geo, oldGeo, newTets, Ynew, Set)
%ADDTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

if ~exist('Ynew', 'var')
    Ynew = [];
end

for newTet = newTets'
    if any(~ismember(newTet, Geo.XgID))
        for numNode = newTet'
            if ~any(ismember(newTet, Geo.XgID)) && ismember(sort(newTet)', sort(Geo.Cells(numNode).T, 2), 'rows')
                Geo.Cells(numNode).Y(ismember(sort(Geo.Cells(numNode).T, 2), sort(newTet)', 'rows'), :) = [];
                Geo.Cells(numNode).T(ismember(sort(Geo.Cells(numNode).T, 2), sort(newTet)', 'rows'), :) = [];
            else
                if isempty(Geo.Cells(numNode).T) || ~ismember(sort(newTet)', sort(Geo.Cells(numNode).T, 2), 'rows')
                    DT = delaunayTriangulation(vertcat(Geo.Cells(newTet).X));
                    if ~any(ismember(newTet, Geo.XgID))
                        Geo.Cells(numNode).T(end+1, :) = newTet;
                    else
                        Geo.Cells(numNode).T(end+1, :) = newTet(DT.ConnectivityList);
                    end
                    if ~isempty(Geo.Cells(numNode).AliveStatus) && exist('Set', 'var')
                        if ~isempty(Ynew)
                            Geo.Cells(numNode).Y(end+1, :) = Ynew(ismember(newTets, newTet', 'rows'), :);
                        else
                            Geo.Cells(numNode).Y(end+1, :) = RecalculateYsFromPrevious(oldGeo, newTet', numNode, Set);
                        end

                        Geo.numY = Geo.numY + 1;
                    end
                end
            end
        end
    end
end
end

