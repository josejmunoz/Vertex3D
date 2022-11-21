function IsConsistent = CheckTris(Geo)
    IsConsistent = true;
    for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        for f = 1:length(Geo.Cells(c).Faces)
            if isempty(vertcat(Geo.Cells(c).Faces(f).Tris.Edge))
                IsConsistent = false;
            end
        end
    end
end