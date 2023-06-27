function [g, K] = initializeKg(Geo, Set)
%%initializeKg 
% 
	dimg=(Geo.numY+Geo.numF+Geo.nCells)*3; 
	g = sparse(dimg, 1);
	K = sparse(dimg, dimg);
end