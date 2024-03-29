function [area, trisArea] = ComputeFaceArea(Tris, Y, FaceCentre)
	area = 0;
    trisArea = cell(size(Tris, 1),1);
	for t = 1:size(Tris, 1)
		Tri = Tris(t,:);
        Y3 = FaceCentre;
		YTri = [Y(Tri,:); Y3];
		T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
        trisArea{t} = T;
		area = area + T;
	end
end