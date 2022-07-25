function [idNeighbours, numNeighbours] = getVertexNeighbours(Geo, idVertex, idCell)
%GETVERTEXNEIGHBOURS Get number and ID of neighbours of a vertex
%   Two tets are neighbours if they share a face.
%   Vertex and Tetrahedra share the same ID

allTets = unique(vertcat(Geo.Cells.T), 'rows');

idNeighbours = find(sum(ismember(allTets, Geo.Cells(idCell).T(idVertex, :)), 2) > 2 & sum(ismember(allTets, Geo.Cells(idCell).T(idVertex, :)), 2) < 4);
numNeighbours = sum(idNeighbours);

end

