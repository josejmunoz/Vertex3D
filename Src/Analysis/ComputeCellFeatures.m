function [features] = ComputeCellFeatures(Cells)
%COMPUTECELLFEATURES Summary of this function goes here
%   Detailed explanation goes here

    for cell = Cells
        if isempty(cell.Area)
            continue
        end
        
        % Compute different measurements from the CELLS
        ComputeCellArea(cell);
        ComputeCellVolume(cell);
        ComputeCellHeight(cell);
        ComputeCellArea(cell, 'Top');
        ComputeCellArea(cell, 'Bottom');
        ComputeCellArea(cell, 'Cell-Cell');
        ComputeCellNeighbours(cell);
        ComputeCellNeighbours(cell, 'Top');
        ComputeCellNeighbours(cell, 'Bottom');
        ComputeCellTilting(cell)
        
        %TODO: Other cell measurements
        %ComputeCellCircularity
        
        % Compute different measurements from the WOUND
        
    end
end

