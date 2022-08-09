% TODO FIXME, should this just be an update to the Geo struct?
function [nrgs]=ComputeTriEnergy(Face, Ys, Set)
    nrgs = zeros(0,1);
	for t = 1:length(Face.Tris)
        y1 = Ys(Face.Tris(t).Edge(1),:);
        y2 = Ys(Face.Tris(t).Edge(2),:);
        y3 = Face.Centre;
        ys(1, :) = {y1, y2, y3};
        ys(2, :) = {y2, y3, y1};
        ys(3, :) = {y3, y1, y2};
        for numY = 1:size(ys, 1)
            y1 = ys{numY, 1}';
            y2 = ys{numY, 2}';
            y3 = ys{numY, 3}';
            [~, cos_angle] = ComputeEdgesAngle(y1, y2, y3);
            cos_angles(numY) = cos_angle;
        end
        nrg = Set.lambdaB/2 * sum((cos_angles - 1/2).^2);
        nrgs(end+1) = nrg;
	end
end
