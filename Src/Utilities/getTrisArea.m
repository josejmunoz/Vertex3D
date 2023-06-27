function faces = getTrisArea(Geo)
	faces = zeros(0,1);
	for c = [Geo.Cells(~cellfun(@isempty, {Geo.Cells.AliveStatus})).ID]
    	for f = 1:length(Geo.Cells(c).Faces)
			trisarea = Geo.Cells(c).Faces(f).TrisArea;
        	faces(end+1:end+length(trisarea),:) = trisarea';
    	end
	end
end