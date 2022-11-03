function [Geo] = AddTetrahedra(Geo, newTets, Ynew, Set)
%ADDTETRAHEDRA Summary of this function goes here
%   Detailed explanation goes here

if ~exist('Ynew', 'var')
    Ynew = []; 
end

for newTet = newTets'
    for numNode = newTet'
        if isempty(Geo.Cells(numNode).T) || all(~ismember(sort(Geo.Cells(numNode).T, 2), sort(newTet', 2), 'rows'))
            DT = delaunayTriangulation(vertcat(Geo.Cells(newTet).X));
            Geo.Cells(numNode).T(end+1, :) = newTet(DT.ConnectivityList);
            if ~isempty(Geo.Cells(numNode).AliveStatus) && exist('Set', 'var')
                if ~isempty(Ynew)
                    Geo.Cells(numNode).Y(end+1, :) = Ynew(ismember(newTets, newTet', 'rows'), :);
                else
                    Geo.Cells(numNode).Y(end+1, :) = ComputeY(Geo, newTet, Geo.Cells(numNode).X, Set);
                end
                
                Geo.numY = Geo.numY + 1;
            end
        end
    end   
end
end

