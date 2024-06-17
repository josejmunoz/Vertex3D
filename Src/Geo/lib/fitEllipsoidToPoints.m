function [a, b, c, paramsOptimized] = fitEllipsoidToPoints(P)
% fitEllipsoidToPoints Fits an ellipsoid to a given set of 3D points.
% The ellipsoid is centered at the origin.
%
% Inputs:
%   P - A 3xN matrix containing the 3D points. Each column represents a point.
%
% Outputs:
%   a, b, c - The semi-axes lengths of the fitted ellipsoid.

    % Extract coordinates from the input matrix
    x = P(:, 1);
    y = P(:, 2);
    z = P(:, 3);

    % Define the objective function for ellipsoid fitting
    function error = ellipsoidError(params)
        % Extract semi-axis lengths from the input parameter
        a = params(1);
        b = params(2);
        c = params(3);

        % Calculate the sum of squared distances from the ellipsoid surface
        distances = (x.^2 / a^2) + (y.^2 / b^2) + (z.^2 / c^2) - 1;
        error = sum(distances.^2);
    end

    % Initial guess for the semi-axis lengths
    initialGuess = [std(x), std(y), std(z)];

    % Set up optimization options
    options = optimoptions('fminunc', 'Algorithm', 'quasi-newton');

    % Perform optimization to find the best-fitting ellipsoid parameters
    paramsOptimized = fminunc(@ellipsoidError, initialGuess, options);

    % Extract optimized parameters
    a = paramsOptimized(1)/max(paramsOptimized);
    b = paramsOptimized(2)/max(paramsOptimized);
    c = paramsOptimized(3)/max(paramsOptimized);
end
