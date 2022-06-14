function [newY] = ComputeY(x, cellCentre, threeGhostNodes)
%COMPUTEY Summary of this function goes here
%   Detailed explanation goes here
    % Condition for the case where 3 nodes are ghost nodes,
    % i.e. external vertex
    if threeGhostNodes
        Center=(sum(x,1))/4;
        vc=Center-cellCentre;
        dir=vc/norm(vc);
        offset=Set.f*dir;
        newY = cellCentre+offset;
    else
        newY = mean(x, [], 2);
    end
end

