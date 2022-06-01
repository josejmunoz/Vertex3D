function [edgeLength] = ComputeEdgeLength(edge, Y)
%COMPUTEEDGELENGTH Summary of this function goes here
%   Detailed explanation goes here
    edgeLength = norm(Y(edge(1), :) - Y(edge(2), :));
end

