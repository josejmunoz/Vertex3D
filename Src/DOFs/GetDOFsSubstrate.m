function [Dofs]=GetDOFsSubstrate(Geo, Set)
    dim = 3;
    gconstrained = zeros((Geo.numY+Geo.numF+Geo.nCells)*3, 1);
    gprescribed  = zeros((Geo.numY+Geo.numF+Geo.nCells)*3, 1);

    for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
        Y     = Geo.Cells(c).Y;
        gIDsY = Geo.Cells(c).globalIds;
        for f = 1:length(Geo.Cells(c).Faces)
            Face = Geo.Cells(c).Faces(f);
            if Face.Centre(3) <= Set.SubstrateZ
            	gconstrained(dim*(Face.globalIds-1)+3) = 1;
            end
        end
        fixY = Y(:,3) <= Set.SubstrateZ;
        for ff = 1:length(find(fixY))
            idx = find(fixY);
            idx = idx(ff);
            gconstrained(dim*(gIDsY(idx)-1)+1:dim*gIDsY(idx)) = 1;
        end
  
    end
    Dofs.Free = find(gconstrained==0 & gprescribed==0);
    Dofs.Fix  = [find(gconstrained); find(gprescribed)];
    Dofs.FixP = find(gprescribed);
    Dofs.FixC = find(gconstrained);
end