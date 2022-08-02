function [] = InitializeGeometry_3DVoronoi()
%INITIALIZEGEOMETRY_3DVORONOI Summary of this function goes here
%   Detailed explanation goes here

nSeeds = 100;
distorsion = 1;
x = rand(nSeeds, 1);
y = rand(nSeeds, 1);

seedsXY = horzcat(x,y);
DT = delaunayTriangulation(x, y);

% figure,
% IC = incenter(DT);
% triplot(DT)
% hold on
% plot(IC(:,1),IC(:,2),'*r')
% hold off

seedsXY_topoChanged = horzcat(x, y + distorsion);
seedsXY_topoChanged(:, 2) = (seedsXY_topoChanged(:, 2) - min(seedsXY_topoChanged(:, 2))) / max(seedsXY_topoChanged(:, 2));

DT_topoChanged = delaunayTriangulation(seedsXY_topoChanged(:, 1), seedsXY_topoChanged(:, 2));
[V, R] = voronoiDiagram(DT_topoChanged);

% figure,
% IC = incenter(DT_topoChanged);
% triplot(DT_topoChanged)
% hold on
% plot(IC(:,1),IC(:,2),'*r')



end

