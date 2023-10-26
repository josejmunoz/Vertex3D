function [features] = ComputeCellFeatures(cell)
%COMPUTECELLFEATURES Summary of this function goes here
%   Detailed explanation goes here
    features = struct();
    % Compute different measurements from the CELLS
    features.Area = ComputeCellArea(cell);
    features.Vol = ComputeCellVolume(cell);
    features.Height = ComputeCellHeight(cell);
    features.Area_Top = ComputeCellArea(cell, 'Top');
    features.Area_Bottom = ComputeCellArea(cell, 'Bottom');
    features.Area_CellCell = ComputeCellArea(cell, 'CellCell');
    features.Neighbours = length(ComputeCellNeighbours(cell));
    features.Neighbours_Top = length(ComputeCellNeighbours(cell, 'Top'));
    features.Neighbours_Bottom = length(ComputeCellNeighbours(cell, 'Bottom'));
    features.Tilting = mean(ComputeCellTilting(cell));
    features.Sphericity = ComputeCellSphericity(cell);
    features.Perimeter_Top = ComputePerimeterByLocation(cell, 'Top');
    features.Perimeter_Bottom = ComputePerimeterByLocation(cell, 'Bottom');
    features.Circularity_Top = Compute2DCircularity(features.Area_Top, features.Perimeter_Top);
    features.Circularity_Bottom = Compute2DCircularity(features.Area_Bottom, features.Perimeter_Bottom);
    [features.Elongation_width, features.Elongation_height, features.Elongation_depth] = ComputeCellElongation(cell);
    [features.Elongation_width_Top, features.Elongation_height_Top, features.Elongation_depth_Top] = ComputeCellElongation(cell, "Top");
    [features.Elongation_width_Bottom, features.Elongation_height_Bottom, features.Elongation_depth_Bottom] = ComputeCellElongation(cell, "Bottom");
    features.AspectRatio2D_Top = Compute2DCellAspectRatio(cell, "Top");
    features.AspectRatio2D_Bottom = Compute2DCellAspectRatio(cell, "Bottom");

    %TODO: Other cell measurements
    %Moments of Inertia: Compute moments of inertia to analyze the object's mass distribution.
    %Descriptors than average non-descript measures like mean length this, maxlenght that, can be obtain by using Chord Length Distributions or Spatial Correlations
    
end

