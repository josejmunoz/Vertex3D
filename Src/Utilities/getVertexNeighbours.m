function [idNeighbours, numNeighbours, tetsNeighbours] = getVertexNeighbours(Geo, idVertex, idCell)
%GETVERTEXNEIGHBOURS Get number and ID of neighbours of a vertex
%   Two tets are neighbours if they share a face.
%   Vertex and Tetrahedra share the same ID

allTets = vertcat(Geo.Cells.T);

idNeighbours = find(sum(ismember(allTets, Geo.Cells(idCell).T(idVertex, :)), 2) == 3);
%Get unique values
[~, uniqueIDs] = unique(sort(allTets(idNeighbours, :), 2), 'rows');
idNeighbours = idNeighbours(uniqueIDs);
numNeighbours = length(idNeighbours);
tetsNeighbours = allTets(idNeighbours, :);

end

