function IsConsistent = CheckTris(Geo)
    IsConsistent = true;
    for c = 1:Geo.nCells
        for f = 1:length(Geo.Cells(c).Faces)
            if isempty(Geo.Cells(c).Faces(f).Tris)
                IsConsistent = false;
            end
        end
    end
end