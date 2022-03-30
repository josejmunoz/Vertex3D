function [Dofs, Geo] = GetRemodelDOFs(Tnew, Dofs, Geo)
	remodelDofs = zeros(0,1);
	for ccc = 1:Geo.nCells
		news = find(sum(ismember(Tnew,ccc)==1,2));
		remodelDofs(end+1:end+length(news),:) = Geo.Cells(ccc).globalIds(end-length(news)+1:end,:);
		for jj = 1:length(Geo.Cells(ccc).Faces)
			Face_r = Geo.Cells(ccc).Faces(jj);
			% TODO FIXME this seems not good...
			FaceTets = Geo.Cells(ccc).T(unique(Face_r.Tris),:);
			if all(ismember(FaceTets,Tnew))
				remodelDofs(end+1,:) = Face_r.globalIds;
			end

		end
	end
	Dofs.Remodel = unique(remodelDofs, 'rows');
	Dofs.Remodel = 3.*(kron(Dofs.Remodel',[1 1 1])-1)+kron(ones(1,length(Dofs.Remodel')),[1 2 3]);
	Geo.AssemblegIds  = unique(unique(remodelDofs, 'rows'));
end