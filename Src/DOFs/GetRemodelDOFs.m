function [Dofs, Geo] = GetRemodelDOFs(Tnew, Dofs, Geo)
remodelDofs = zeros(0,1);
for numCell = 1:Geo.nCells
    news = ismember(sort(Geo.Cells(numCell).T, 2), sort(Tnew, 2), 'rows');
    remodelDofs(end+1:end+sum(news),:) = Geo.Cells(numCell).globalIds(news);
    for jj = 1:length(Geo.Cells(numCell).Faces)
        Face_r = Geo.Cells(numCell).Faces(jj);
        FaceTets = Geo.Cells(numCell).T(unique([Face_r.Tris.Edge]),:);
        if any(ismember(sort(FaceTets, 2),sort(Tnew, 2), 'rows'))
            remodelDofs(end+1,:) = Face_r.globalIds;
        end
    end
end
Dofs.Remodel = unique(remodelDofs, 'rows');
Dofs.Remodel = 3.*(kron(Dofs.Remodel',[1 1 1])-1)+kron(ones(1,length(Dofs.Remodel')),[1 2 3]);
Geo.AssemblegIds  = unique(unique(remodelDofs, 'rows'));
end