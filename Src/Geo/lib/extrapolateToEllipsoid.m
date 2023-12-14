function [points_ellipsoid] = extrapolateToEllipsoid(points, a, b, c)
%EXTRAPOLATETOELLIPSOID Summary of this function goes here
%   Detailed explanation goes here
% Assuming you have a set of points on the sphere stored in a matrix 'points'
% Each row of 'points' contains the (x, y, z) coordinates of a point

    % Define the number of points
    points_ellipsoid = [points(:, 1) * a, points(:, 2) * b, points(:, 3) * c];
    
%     % Plot the points on the original sphere-like surface and the extrapolated ellipsoid
%     figure;
%     scatter3(points(:, 1), points(:, 2), points(:, 3), 10, 'b', 'filled'); % Points on the sphere
%     hold on;
%     scatter3(points_ellipsoid(:, 1), points_ellipsoid(:, 2), points_ellipsoid(:, 3), 10, 'r', 'filled'); % Extrapolated ellipsoid points
end