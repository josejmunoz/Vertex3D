function [nrgs]=ComputeTriAREnergy(Face, Ys, Set)
    nrgs = zeros(0,1);
	for t = 1:length(Face.Tris)
        y1 = Ys(Face.Tris(t).Edge(1),:);
        y2 = Ys(Face.Tris(t).Edge(2),:);
        y3 = Face.Centre;
        ys(1, :) = {y1, y2, y3};
        ys(2, :) = {y2, y3, y1};
        ys(3, :) = {y3, y1, y2};
        w_t = zeros(3, 1);
        for numY = 1:size(ys, 1)
            y1 = ys{numY, 1}';
            y2 = ys{numY, 2}';
            y3 = ys{numY, 3}';

            v_y1 = y2 - y1;
            v_y2 = y3 - y1;

            w_t(numY) = norm(v_y1)^2 - norm(v_y2)^2;
        end
        nrg = Set.lambdaR/2 * sum(w_t.^2) * 1/(Set.lmin0^4);
        nrgs(end+1) = nrg;
	end
end
