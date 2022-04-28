function Geo = ReplaceYs(targetTets, Tnew, Ynew, Geo)

	targetNodes = unique(targetTets);
	for n_i = 1:length(targetNodes)
		tNode = targetNodes(n_i);
		CellJ = Geo.Cells(tNode);
		hits = find(sum(ismember(CellJ.T,targetTets),2)==4);
		Geo.Cells(tNode).T(hits,:) = [];

		news = find(sum(ismember(Tnew,tNode)==1,2));
		Geo.Cells(tNode).T(end+1:end+length(news),:) = Tnew(news,:);
		if ~ismember(tNode, Geo.XgID)
			Geo.Cells(tNode).Y(hits,:) = [];
			Geo.Cells(tNode).Y(end+1:end+length(news),:) = Ynew(news,:);
		end
	end
end