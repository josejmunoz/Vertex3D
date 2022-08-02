function [trianglesConnectivity, neighboursNetwork, cellEdges_Boundary, verticesOfCell_pos] = Build3DVoronoiTopo(seedsXY)
%BUILD3DVORONOITOPO Summary of this function goes here
%   Detailed explanation goes here
DT = delaunayTriangulation(seedsXY(:, 1), seedsXY(:, 2));
triangleNeighbours = neighbors(DT); %% Connected to DT.ConnectivityList
%cellEdges = vertexAttachments(DT);
neighboursNetwork = edges(DT); %% OK
trianglesConnectivity = DT.ConnectivityList; %% OK
verticesOfCell_pos = circumcenter(DT); %% OK

cellEdges_Boundary = cell(1, size(seedsXY, 1));
for numCell = 1:size(seedsXY, 1)
    verticesIndices = find(any(ismember(DT.ConnectivityList, numCell), 2));
    
    if length(verticesIndices) > 2
        cellEdgesOrdered = verticesIndices(boundary(verticesOfCell_pos(verticesIndices, :)));
        cellEdgesOrdered(:, 2) = cellEdgesOrdered([2:end 1]);
        cellEdges_Boundary(numCell) = {cellEdgesOrdered(1:end-1, :)};
    end
end

% figure,
% IC = incenter(DT);
% triplot(DT)
% hold on
% plot(IC(:,1),IC(:,2),'*r')
% hold off
end

