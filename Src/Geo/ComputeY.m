function [newY] = ComputeY(Geo, T, cellCentre, Set)
%COMPUTEY Summary of this function goes here
%   Detailed explanation goes here
    % Condition for the case where 3 nodes are ghost nodes,
    % i.e. external vertex
    
    x = vertcat(Geo.Cells(T).X);
    newY = mean(x);
    
    if length([Geo.Cells(T).AliveStatus]) == 1 && contains(Set.InputGeo, 'Bubbles')
        vc=newY-cellCentre;
        dir=vc/norm(vc);
        offset=Set.f*dir;
        newY = cellCentre+offset;
    end
    
    if ~contains(Set.InputGeo, 'Bubbles')
        if sum(ismember(T, Geo.XgTop)) > 0
            newY(3) = newY(3) / (sum(ismember(T, Geo.XgTop))/2);
        elseif sum(ismember(T, Geo.XgBottom)) > 0
            newY(3) = newY(3) / (sum(ismember(T, Geo.XgBottom))/2);
        end
    end

    % function moveVertexToSurfaceOfCells
    if Set.TryingFlips == 0
        if sum(arrayfun(@(x) size(x.T, 1) == size(x.Y, 1), Geo.Cells(T))) == sum(~cellfun(@isempty, {Geo.Cells(T).AliveStatus}))
            newY = moveVertexToSurfaceOfCell(Geo, T, newY);
        end
    end
end

