function [newVertOrder] = BoundaryOfCell(verticesOfCell, neighbours)
%BOUNDARYOFCELL Summary of this function goes here
% Here we consider 3 methods to connect the vertices, and we choose the
% method with more area into the polyshape.
% By Pedro J. Gomez Galvez, modified

if exist('neighbours', 'var')
    try
        initialNeighbours = neighbours;
        neighboursOrder = neighbours(1, :);
        firstNeighbour = neighbours(1, 1);
        nextNeighbour = neighbours(1, 2);
        nextNeighbourPrev = nextNeighbour;
        neighbours(1, :) = [];
        while isempty(neighbours) == 0
            matchNextVertex = any(ismember(neighbours, nextNeighbour), 2);

            neighboursOrder(end+1, :) = neighbours(matchNextVertex, :);

            nextNeighbour = neighbours(matchNextVertex, :);
            nextNeighbour(nextNeighbour == nextNeighbourPrev) = [];
            neighbours(matchNextVertex, :) = [];

            nextNeighbourPrev = nextNeighbour;
        end

        [~, vertOrder] = ismember(neighboursOrder, initialNeighbours, 'rows');
        
        newVertOrder = horzcat(vertOrder, vertcat(vertOrder(2:end), vertOrder(1)));
        
        return;
    catch ex
        newVertOrder = [];
    end
end

imaginaryCentroidMeanVert = mean(verticesOfCell);
vectorForAngMean = bsxfun(@minus, verticesOfCell, imaginaryCentroidMeanVert );
thMean = atan2(vectorForAngMean(:,2),vectorForAngMean(:,1));
[~, vertOrder] = sort(thMean);

newVertOrder = horzcat(vertOrder, vertcat(vertOrder(2:end), vertOrder(1)));
end