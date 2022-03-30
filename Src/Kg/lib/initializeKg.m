function [g, K] = initializeKg(Geo, Set)
		dimg=(Geo.numY+Geo.numF+Geo.nCells)*3; 

	g = zeros(dimg, 1);
	K = zeros(dimg, dimg);
end