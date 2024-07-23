function [tilting] = ComputeEdgeTilting(edge, Y)
%COMPUTEEDGETILTING Summary of this function goes here
%   Detailed explanation goes here
    v1 = Y(edge.Edge(1), :) - Y(edge.Edge(2), :); % realEdge
    
    if edge.Location == 1
        fixedVertex = [Y(edge.Edge(1), 1:2), Y(edge.Edge(2), 3)];
    else
        fixedVertex = [Y(edge.Edge(2), 1:2), Y(edge.Edge(1), 3)];
    end
    %TODO: CHECK IF THIS IS CORRECT
    %TODO: Improve perpendicular edge for curve tissues
    v2 = Y(edge.Edge(1), :) - fixedVertex;% Perpendicular edge
    tilting = atan2(norm(cross(v1,v2)),dot(v1,v2)) * 100;
end

