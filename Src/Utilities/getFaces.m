function faces = getFaces(Geo)
	faces = zeros(0,3);
	for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
    	for f = 1:length(Geo.Cells(c).Faces)
        	faces(end+1,:) = Geo.Cells(c).Faces(f).Centre;
    	end
	end
	faces = unique(faces,'rows');
end