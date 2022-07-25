function [idNeighbours, numNeighbours] = getVertexNeighbours(Geo, idVertex, idCell)
%GETVERTEXNEIGHBOURS Get number and ID of neighbours of a vertex
%   Two tets are neighbours if they share a face.
%   Vertex and Tetrahedra share the same ID

allTets = unique(vertcat(Geo.Cells.T), 'rows');

idNeighbours = sum(ismember(allTets, Geo.Cells(idCell).T(idVertex, :)), 2) > 2;
idNeighbours(trisToChange.Edge(1)) = 0;
numNeighbours = sum(idNeighbours);


end

