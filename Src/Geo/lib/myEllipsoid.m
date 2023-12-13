function [x_ellipsoid, y_ellipsoid, z_ellipsoid] = myEllipsoid(a, b, c, num_points)
%MYELLIPSOID Summary of this function goes here
%   
%     a = Semi-axis length along x-axis
%     b = Semi-axis length along y-axis
%     c = Semi-axis length along z-axis

    % Define the parameters of the ellipsoid
    x0 = 0; % Center x-coordinate
    y0 = 0; % Center y-coordinate
    z0 = 0; % Center z-coordinate

    % Normalise axis
    a = a/max([a,b,c]);
    b = b/max([a,b,c]);
    c = c/max([a,b,c]);
    
    % Generate random points within a unit sphere
    theta = 2 * pi * rand(num_points, 1);
    phi = acos(2 * rand(num_points, 1) - 1);
    r = 1; % Unit sphere radius
    
    % Convert spherical coordinates to Cartesian coordinates
    x_sphere = r * sin(phi) .* cos(theta);
    y_sphere = r * sin(phi) .* sin(theta);
    z_sphere = r * cos(phi);
    
    % Scale the points to fit the ellipsoid
    x_ellipsoid = x0 + a * x_sphere;
    y_ellipsoid = y0 + b * y_sphere;
    z_ellipsoid = z0 + c * z_sphere;
end

