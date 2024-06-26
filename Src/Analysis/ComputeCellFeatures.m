function [features] = ComputeCellFeatures(cell)
%COMPUTECELLFEATURES Summary of this function goes here
%   Detailed explanation goes here
    features = struct();
    % Compute different measurements from the CELLS
    features.Area = ComputeCellArea(cell);
    features.Vol = ComputeCellVolume(cell);
    features.Height = ComputeCellHeight(cell);
    features.Area_Top = ComputeCellArea(cell, 1);
    features.Area_Bottom = ComputeCellArea(cell, 3);
    features.Area_CellCell = ComputeCellArea(cell, 2);
    features.Neighbours = length(ComputeCellNeighbours(cell));
    features.Neighbours_Top = length(ComputeCellNeighbours(cell, 1));
    features.Neighbours_Bottom = length(ComputeCellNeighbours(cell, 3));
    features.Tilting = mean(ComputeCellTilting(cell));
    features.Sphericity = ComputeCellSphericity(cell);
    features.Perimeter_Top = ComputePerimeterByLocation(cell, 1);
    features.Perimeter_Bottom = ComputePerimeterByLocation(cell, 3);
    features.Circularity_Top = Compute2DCircularity(features.Area_Top, features.Perimeter_Top);
    features.Circularity_Bottom = Compute2DCircularity(features.Area_Bottom, features.Perimeter_Bottom);
    [features.Elongation_width, features.Elongation_height, features.Elongation_depth] = ComputeCellElongation(cell);
    [features.Elongation_width_Top, features.Elongation_height_Top, features.Elongation_depth_Top] = ComputeCellElongation(cell, 1);
    [features.Elongation_width_Bottom, features.Elongation_height_Bottom, features.Elongation_depth_Bottom] = ComputeCellElongation(cell, 3);
    features.AspectRatio2D_Top = Compute2DCellAspectRatio(cell, 1);
    features.AspectRatio2D_Bottom = Compute2DCellAspectRatio(cell, 3);

    %TODO: Other cell measurements
    %Moments of Inertia: Compute moments of inertia to analyze the object's mass distribution.
    %Descriptors than average non-descript measures like mean length this, maxlenght that, can be obtain by using Chord Length Distributions or Spatial Correlations
    
end

