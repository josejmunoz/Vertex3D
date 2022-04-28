function [Geo]=BuildCells(Geo, Set)
	% Iterate over cells
	for c = 1:length(Geo.Cells)
		Tets = Geo.Cells(c).T;
		for f = 1:length(Tets)
			for v = 1:nverts
				
			end
		end
	end
end