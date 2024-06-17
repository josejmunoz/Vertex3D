function [width, height, depth] = ComputeCellElongation(cell, location)
%COMPUTECELLELONGATION Summary of this function goes here
%   Detailed explanation goes here
% Getting x and y data

try
    if ~exist('location', 'var')
        x = cell.Y(:, 1);
        y = cell.Y(:, 2);
        z = cell.Y(:, 3);
    else
        vertices = [];
        for face = cell.Faces
            for tri = face.Tris
                if face.InterfaceType == location
                    vertices = vertcat(vertices, vertcat(cell.Y(tri.Edge, :)));
                end
            end
        end
        x = vertices(:, 1);
        y = vertices(:, 2);
        z = vertices(:, 3);
    end
    % Getting width and height
    width = max(x) - min(x);
    height = max(y) - min(y);
    depth = max(z) - min(z);
catch
    disp('Possible error on ComputeCellElongation')
    width = -1;
    height = -1;
    depth = -1;
end
end

