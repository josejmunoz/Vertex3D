function [newY] = ComputeY(x, cellCentre, lengthOfAlivedCells, Set)
%COMPUTEY Summary of this function goes here
%   Detailed explanation goes here
    % Condition for the case where 3 nodes are ghost nodes,
    % i.e. external vertex
    newY = mean(x);
    if lengthOfAlivedCells == 1
        vc=newY-cellCentre;
        dir=vc/norm(vc);
        offset=Set.f*dir;
        newY = cellCentre+offset;
    end
end

