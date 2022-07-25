function [Geo_n, Geo, Dofs, Set, newYgIds] = Flip13(Geo_0, Geo_n, Geo, Dofs, Set, newYgIds)
%FLIP13 Summary of this function goes here
%   Detailed explanation goes here

for c = 1:Geo.nCells
    
    f = 0;
    allTris = [Geo.Cells(c).Faces.Tris];
    avgArea = mean([allTris.Area]);
    stdArea = std([allTris.Area]);
    
    %CARE: Number of faces change within this loop, so it should be a while
    while f < length(Geo.Cells(c).Faces)
        f = f + 1;
        
        Ys = Geo.Cells(c).Y;
        Ts = Geo.Cells(c).T;
        
        Face = Geo.Cells(c).Faces(f);
        faceAreas = [Face.Tris.Area];
        [maxTriArea, idMaxTriArea]= max(faceAreas);
        
        if maxTriArea < avgArea + stdArea*2 || ismember(Face.globalIds, newYgIds)
            continue
        end
        
        Geo_backup = Geo; Geo_n_backup = Geo_n;
        
        trisToChange = Face.Tris(idMaxTriArea);

        [idNeighbours_1, numNeighbours_1] = getVertexNeighbours(Geo, trisToChange.Edge(1), c);
        [idNeighbours_2, numNeighbours_2] = getVertexNeighbours(Geo, trisToChange.Edge(2), c);
        
        if numNeighbours_1 == 3 || numNeighbours_2 == 3
            idNeighbours_1
        end
    end
end
end

