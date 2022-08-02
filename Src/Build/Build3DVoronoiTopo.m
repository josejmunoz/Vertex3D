function [trianglesConnectivity, neighboursNetwork, cellEdges, verticesOfCell_pos] = Build3DVoronoiTopo(seedsXY)
%BUILD3DVORONOITOPO Summary of this function goes here
%   Detailed explanation goes here
DT = delaunayTriangulation(seedsXY(:, 1), seedsXY(:, 2));
triangleNeighbours = neighbors(DT); %% Connected to DT.ConnectivityList
cellEdges = vertexAttachments(DT);
neighboursNetwork = edges(DT);
trianglesConnectivity = DT.ConnectivityList;

[verticesOfCell_pos] = voronoiDiagram(DT);
verticesOfCell_pos = verticesOfCell_pos(2:end, :);

% figure,
% IC = incenter(DT);
% triplot(DT)
% hold on
% plot(IC(:,1),IC(:,2),'*r')
% hold off
end

