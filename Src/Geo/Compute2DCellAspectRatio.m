function [] = Compute2DCellAspectRatio(cell, location)
%COMPUTE2DCELLASPECTRATIO Summary of this function goes here
%   Detailed explanation goes here
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

ellipse_t = fit_ellipse( x,y );
end

