function [aspectRatio] = Compute2DCellAspectRatio(cell, location)
%COMPUTE2DCELLASPECTRATIO Summary of this function goes here
%   Detailed explanation goes here

try
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
    if ~isempty(ellipse_t)
        aspectRatio = ellipse_t.long_axis/ellipse_t.short_axis;
    else
        aspectRatio = 0;
    end
catch
    disp('Possible error on ComputeCell2DAspectRatio')
    aspectRatio = -1;
end
end

