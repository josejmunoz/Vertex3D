function [Dofs, Geo] = GetRemodelDOFs(Tnew, Dofs, Geo)
remodelDofs = zeros(0,1);
aliveCells = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID];
idTnew = unique(Tnew);
idTnew_cells = setdiff(idTnew, Geo.XgID);
idTnew_cells = idTnew_cells([Geo.Cells(idTnew_cells).AliveStatus] == 1);
for numCell = idTnew_cells'
    news = sum(ismember(Geo.Cells(numCell).T, Geo.XgID), 2) > 2;
    news(sum(ismember(Geo.Cells(numCell).T, idTnew_cells), 2) == 2 & sum(ismember(Geo.Cells(numCell).T, Geo.XgID), 2) == 2) = 1;
    news(sum(ismember(Geo.Cells(numCell).T, idTnew_cells), 2) >= 3) = 1;

    % Remove only the tets from the domain it is not changing
    if sum(sum(ismember(Tnew, Geo.XgBottom))) > sum(sum(ismember(Tnew, Geo.XgTop)))
        news(any(ismember(Geo.Cells(numCell).T, Geo.XgTop), 2)) = 0;
    else
        news(any(ismember(Geo.Cells(numCell).T, Geo.XgBottom), 2)) = 0;
    end

    remodelDofs(end+1:end+sum(news),:) = Geo.Cells(numCell).globalIds(news);
    for jj = 1:length(Geo.Cells(numCell).Faces)
        Face_r = Geo.Cells(numCell).Faces(jj);
        FaceTets = Geo.Cells(numCell).T(unique([Face_r.Tris.Edge]),:);
        if any(ismember(sort(FaceTets, 2),sort(Geo.Cells(numCell).T(news, :), 2), 'rows'))
            remodelDofs(end+1,:) = Face_r.globalIds;
        end
    end
end
Dofs.Remodel = unique(remodelDofs, 'rows');
Dofs.Remodel = 3.*(kron(Dofs.Remodel',[1 1 1])-1)+kron(ones(1,length(Dofs.Remodel')),[1 2 3]);
Geo.AssemblegIds  = unique(unique(remodelDofs, 'rows'));
end