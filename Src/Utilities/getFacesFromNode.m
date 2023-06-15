function [faces, facesTris] = getFacesFromNode(Geo, nodes)
%GETFACESFROMNODE Summary of this function goes here
%   Detailed explanation goes here
    faces = {};
    for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        for f = 1:length(Geo.Cells(c).Faces)
            if all(ismember(nodes, Geo.Cells(c).Faces(f).ij))
                faces{end+1} = Geo.Cells(c).Faces(f);
            end
        end
    end
    
    allFaces = [faces{:}];
    facesTris = [allFaces(:).Tris];
end

