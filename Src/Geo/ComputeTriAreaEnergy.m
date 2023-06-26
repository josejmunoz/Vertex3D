function [nrgs]=ComputeTriAreaEnergy(Face, Set)
    nrgs = zeros(0,1);
	for t = 1:length(Face.Tris)
        nrg = exp(Set.lambdaB*(1-Set.Beta*Face.Tris(t).Area/Set.BarrierTri0));
        nrgs(end+1) = nrg;
	end
end
