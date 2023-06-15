function [perimeter, trisPerimeter] = ComputeFacePerimeter(Tris, Y, FaceCentre)
%COMPUTEFACEPERIMETER Summary of this function goes here
%   Detailed explanation goes here
	perimeter = 0;
    trisPerimeter = cell(length(Tris),1);
	for t = 1:length(Tris)
		Tri = Tris(t,:);
        Y3 = FaceCentre;
		YTri = [Y(Tri,:); Y3];
		T = norm(YTri(1, :) - YTri(2, :)) + norm(YTri(2, :) - YTri(3, :)) + norm(YTri(3, :) - YTri(1, :));
        trisPerimeter{t} = T;
		perimeter = perimeter + T;
	end
end

