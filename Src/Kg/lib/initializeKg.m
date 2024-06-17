function [g, K] = initializeKg(Geo, Set)
%%initializeKg 
% 
	dimg=(Geo.numY+Geo.numF+Geo.nCells)*3; 
	g = spalloc(dimg, 1, dimg);
	K = spalloc(dimg, dimg, round(dimg*(dimg/10)));
	%K = gpuArray(zeros(dimg, dimg));
end