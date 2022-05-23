function [area, trisArea] = ComputeFaceArea(Tris, Y, FaceCentre)
	area = 0;
    trisArea = cell(length(Tris),1);
	for t = 1:length(Tris)
		Tri = Tris(t,:);

        if length(Tris)==3
            Y3 = Y(Tris(t+1,2),:);
        else
            Y3 = FaceCentre;
        end
		YTri = [Y(Tri,:); Y3];
		T=(1/2)*norm(cross(YTri(2,:)-YTri(1,:),YTri(1,:)-YTri(3,:)));
        trisArea{t} = T;
		area = area + T;
        if length(Tris)==3
            break
        end
	end
end