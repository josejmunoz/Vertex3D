function [trisArea, trisEdges] = getTrisAreaOfNeighbours(Geo, vertexToAnalyse)
%GETTRISAREAOFNEIGHBOURS Summary of this function goes here
%   vertexToAnalyse has to be globalID

trisArea = [];
trisEdges = [];
for cCell = Geo.Cells
    for face = cCell.Faces
        for tris = face.Tris
            if ismember(vertexToAnalyse, cCell.globalIds(tris.Edge))
                trisArea(end+1) = tris.Area;
                trisEdges(end+1, 1:2) = cCell.globalIds(tris.Edge);
            end
        end
    end
end
end

