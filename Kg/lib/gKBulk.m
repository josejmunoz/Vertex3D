function [g, K] = gKBulk(y1, y2, y3, y4)
%GKBULK Summary of this function goes here
%   INPUT:
%   y1, y2, y3, y4 belongs to a single tetrahedron
%
%   OUTPUT:
%   g: 4x3 matrix (4 vertices x 3 coordinates)
%   K: 12x12 matrix

dim = 3;

g = zeros(4, dim);

K = zeros(numel(g), numel(g));

end

