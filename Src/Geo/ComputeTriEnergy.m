% TODO FIXME, should this just be an update to the Geo struct?
function [nrgs]=ComputeTriEnergy(Face, Ys, Set)
    nrgs = zeros(0,1);
	for t = 1:length(Face.Tris)
        %Tri = Face.Tris(t).Edge;
        %Y3 = Face.Centre;
        %YTri = [Ys(Tri(1),:); Ys(Tri(2),:); Y3];
        %area = (1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
        %nrg  = exp(Set.lambdaB*(1-Set.Beta*area/Set.BarrierTri0));
        [sideLengths] = ComputeTriSideLengths(Face, t, Ys);
        [aspectRatio] = ComputeTriAspectRatio(sideLengths);
        nrg  = exp(Set.lambdaB*6*(1/(1-aspectRatio)));
%         aspectRatio
%         nrg
        nrgs(end+1) = nrg;
	end
end
