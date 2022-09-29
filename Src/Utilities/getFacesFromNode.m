function [faces, facesTris] = getFacesFromNode(Geo, node)
%GETFACESFROMNODE Summary of this function goes here
%   Detailed explanation goes here
    faces = {};
    for c = 1:Geo.nCells
        for f = 1:length(Geo.Cells(c).Faces)
            if any(ismember(Geo.Cells(c).Faces(f).ij, node))
                faces{end+1} = Geo.Cells(c).Faces(f);
            end
        end
    end
    
    allFaces = [faces{:}];
    facesTris = [allFaces(:).Tris];
end

